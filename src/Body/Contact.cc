#include "Body.h"
#include "Simulation/Simulation.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "Diagnostics/Operation_Count.h"
#include "List.h"
#include "Array.h"
#include <math.h>

struct Contact_Particle_Bucket {
  Array<unsigned> Array_A{};
  Array<unsigned> Array_B{};

  unsigned Counter_A = 0;
  unsigned Counter_B = 0;
}; // typdef struct Particle_Bucket {

/* Global (shared) variables for the buckets. */
static Contact_Particle_Bucket * Buckets;        // Array of buckets
static unsigned * Bucket_Indicies_Body_A;        // Which bucket each particle of body A goes to.
static unsigned * Bucket_Indicies_Body_B;        // Which bucket each paritcle of body B goes to.
static double Buffer[6];


void Body::Contact(Body & Body_A, Body & Body_B) {
  /* This function calculates the contact force between the particles in both A
  and those in body B.

  First, we partition the spatial domain. To do this, we first determine the
  maximum and minimum x, y, and z coordinates of the live (damage < 1) particles
  This gives us a box in which the particles of the two bodies live. We
  then divide the x, y, and z dimensions of this box into smaller boxes,
  called cells, each one of which has a side length that is barely greater than
  the contact distance. We then allocate an array of buckets with one bucket per
  cell. We then cycle through the particles of A and B, determining which bucket
  each particle belongs to.

  Once we finish this, we allocate Array_A and array_B in each bucket. We
  populate Array_A in the ith bucket with the ID's of the particles in body A
  that live in that bucket. Once these arrays have been populated, the particles
  have essentially been "sorted" into their buckets.

  We then apply the contact algorithm. The way that these buckets are
  set up, a particle can only contact particles that are in its bucket or a
  bucket which is adjacent to its bucket (27 total buckets to check). This
  greatly reduces the number of particles that we need to check against, thereby
  eliminating needless calculations. */



  //////////////////////////////////////////////////////////////////////////////
  /* Determine how many buckets we need */

  // Get parameters from Body_A, Body_B
  const unsigned Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned Num_Particles_B = Body_B.Get_Num_Particles();

  // Declare some variables.
  double x_max, x_min;
  double y_max, y_min;
  double z_max, z_min;

  // Initialize the variables.
  Vector x = Body_A[0].Get_x();
  x_max = x[0];
  x_min = x_max;
  y_max = x[1];
  y_min = y_max;
  z_max = x[2];
  z_min = z_max;

  /* Determine the maximum and minimum x, y, z coordinate of the particles
  in the two bodies. Each thread only searches through a subset of the
  particles. Once this is done, we will "reduce" them together in the
  buffer to determine the global minimum. */
  #pragma omp for nowait
  for(unsigned i = 0; i < Num_Particles_A; i++) {
    x = Body_A[i].Get_x();

    /* Note: if x[0] > x_max, then we can't also have x[0] < x_min (this relies
    on the fact that x_min <= x_max) */
    if     (x[0] > x_max) { x_max = x[0]; }
    else if(x[0] < x_min) { x_min = x[0]; }

    if     (x[1] > y_max) { y_max = x[1]; }
    else if(x[1] < y_min) { y_min = x[1]; }

    if     (x[2] > z_max) { z_max = x[2]; }
    else if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  #pragma omp for nowait
  for(unsigned i = 0; i < Num_Particles_B; i++) {
    x = Body_B[i].Get_x();

    if     (x[0] > x_max) { x_max = x[0]; }
    else if(x[0] < x_min) { x_min = x[0]; }

    if     (x[1] > y_max) { y_max = x[1]; }
    else if(x[1] < y_min) { y_min = x[1]; }

    if     (x[2] > z_max) { z_max = x[2]; }
    else if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  /* The static global buffer variable will be used to perform the reduce
  operation. We will associate the 6 elements of the Buffer to perform the
  following reductions:
  Buffer[0] = global x_min
  Buffer[1] = global x_max
  Buffer[2] = global y_min
  Buffer[3] = global y_max
  Buffer[4] = global z_min
  Buffer[5] = global z_max */

  /* First, we need to initialize the buffer. Whichever thread gets here first
  does that. We need the implict barrier so that threads don't begin the
  reduction until the buffer has been initialized. */
  #pragma omp single
  {
    Buffer[0] = x_min;
    Buffer[1] = x_max;
    Buffer[2] = y_min;
    Buffer[3] = y_max;
    Buffer[4] = z_min;
    Buffer[5] = z_max;
  } // #pragma omp single

  /* Now, perform the reduction. */
  #pragma omp critical
  {
    if(x_min < Buffer[0]) { Buffer[0] = x_min; }
    if(x_max > Buffer[1]) { Buffer[1] = x_max; }

    if(y_min < Buffer[2]) { Buffer[2] = y_min; }
    if(y_max > Buffer[3]) { Buffer[3] = y_max; }

    if(z_min < Buffer[4]) { Buffer[4] = z_min; }
    if(z_max > Buffer[5]) { Buffer[5] = z_max; }
  } // #pragma omp critical
  #pragma omp barrier

  /* Once the reduction is done, each thread can get the global min/max
  quantities from the buffer. */
  x_min = Buffer[0];
  x_max = Buffer[1];
  y_min = Buffer[2];
  y_max = Buffer[3];
  z_min = Buffer[4];
  z_max = Buffer[5];

  #ifdef CONTACT_MONITOR
    /* Check that every thread got the same min/max x, y, and z values. */
    #pragma omp critical
    {
      printf("Buffer = %p\n", Buffer);
      printf("x_min = %lf\n", x_min);
      printf("x_max = %lf\n", x_max);
      printf("y_min = %lf\n", y_min);
      printf("y_max = %lf\n", y_max);
      printf("z_min = %lf\n", z_min);
      printf("z_max = %lf\n", z_max);
    } // #pragma omp critical
  #endif



  //////////////////////////////////////////////////////////////////////////////
  /* Now, determine the bucket dimensions and allocate the buckets */

  /* We want the dimension (in all three coordinate directions) of the
  cell to be >= the contact distance. By doing this, particles
  can only compe into contact with particles in their bucket or in buckets
  that are adjacent to their bucket. In general, we want the buckets to be as
  small as possible (so that there are as few particles to check for contact as
  possible). Let's focus on the x coordinate. Let Nx denote the number of
  cell in the x direction. We want Nx to be the largest natural number
  such that Contact_Distance <= (x_max - x_min)/Nx. A little though reveals
  that this occurs precisely when Nx = floor((x_max - x_min)/Contact_Distance).
  A similar result holds for the y and z directions. */
  const unsigned Nx = floor((x_max - x_min)/Simulation::Contact_Distance);
  const unsigned Ny = floor((y_max - y_min)/Simulation::Contact_Distance);
  const unsigned Nz = floor((z_max - z_min)/Simulation::Contact_Distance);

  #ifdef OPERATION_COUNT
    // 3 subtractions, 3 divisions in the computations above.
    #pragma omp single nowait
    {
      #pragma omp atomic update
      OP_Count::Subtraction += 3;
      #pragma omp atomic update
      OP_Count::Division += 3;
    } // #pragma omp single nowait
  #endif

  /* Now that we know the number of buckets in each direction, one thread can
  allocate the global arrays */
  #pragma omp single
  {
    Buckets = new Contact_Particle_Bucket[Nx*Ny*Nz];
    Bucket_Indicies_Body_A = new unsigned[Num_Particles_A];
    Bucket_Indicies_Body_B = new unsigned[Num_Particles_B];
  } // #pragma omp single

  #ifdef CONTACT_MONITOR
    /* Check that Nx, Ny, Nz were set up properly and that all threads see the
    same dynamically allocated global arrays */
    #pragma omp critical
    {
      printf("(Nx, Ny, Nz) = (%d, %d, %d)\n", Nx, Ny, Nz);
      printf("Buckets = %p\n", Buckets);
      printf("Bucket_Indicies_Body_A = %p\n", Bucket_Indicies_Body_A);
      printf("Bucket_Indicies_Body_B = %p\n", Bucket_Indicies_Body_B);
    } // #pragma omp critical
  #endif



  //////////////////////////////////////////////////////////////////////////////
  /* Now, determine which bucket each particle belongs to. */

  /* First, we calculate some variables (see the next comment for an
  explanation) */
  const double cell_x_dim = (x_max - x_min)/Nx;
  const double cell_y_dim = (y_max - y_min)/Ny;
  const double cell_z_dim = (z_max - z_min)/Nz;

  #ifdef OPERATION_COUNT
    // 3 subtractions, 3 divisions in the computations above.
    #pragma omp single nowait
    {
      #pragma omp atomic update
      OP_Count::Subtraction += 3;
      #pragma omp atomic update
      OP_Count::Division += 3;
    } // #pragma omp single nowait
  #endif

  #pragma omp for nowait
  for(unsigned i = 0; i < Num_Particles_A; i++) {
    /* First, we need to determine which bucket our particle belongs in. Let's
    focus on the x coordinate. Each bucket has an x-dimension length of
    (x_max - x_min)/Nx, which we call cell_x_dim. The x coordinates of
    the bucket for the ith particle is the number of units of length
    cell_x_dim that fit between x_min and the x coordinate of the particle
    (think about it). A little thought reveals that this is precisely
    floor((Particle[i].x[0] - x_min)/cell_x_dim)
    A similar result holds for y and z. */
    x = Body_A.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/cell_x_dim);
    unsigned ny = floor((x[1] - y_min)/cell_y_dim);
    unsigned nz = floor((x[2] - z_min)/cell_z_dim);

    /* We run into a bit of a problem if a pritlce's x coordinate is equal to
    x_max. In this case (assuming no roundoff error), nx will evaluate to Nx.
    This is problematic, because the bucket x coordinates range from 0 to Nx-1
    (remember, 0 indexing). Really, we want this particle to go into a bucket
    corresonding to the cell with x index Nx-1 (think about it). To remedy this,
    we simply run a check: if nx evaluated to Nx, then correct nx to Nx-1.

    This is quite literaly an edge case. It's possible for particles to lie on
    the boundary of other cells. In those cases, however, it's fine to assign
    the particle to either cell (think about it). Thus, we only really need to
    run the check when the code tries to assign a particle to a cell with too
    big an index.

    A similar argument holds for ny and nz. */
    if(nx == Nx) { nx = Nx - 1; }
    if(ny == Ny) { ny = Ny - 1; }
    if(nz == Nz) { nz = Nz - 1; }

    /* Determine which bucket this particle goes into. Also, increment the
    number of particles of body A that go in that bucket. Note that we need
    the atomic declaration because it is theoretically possible for two
    threads to increment the same bucket element at the same time. */
    Bucket_Indicies_Body_A[i] = nx + ny*Nx + nz*Nx*Ny;

    /* Update that bucket's A counter */
    #pragma omp atomic
    Buckets[nx + ny*Nx + nz*Nx*Ny].Counter_A++;
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  #pragma omp for
  for(unsigned i = 0; i < Num_Particles_B; i++) {
    /* Calculate bucket indicies */
    x = Body_B.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/cell_x_dim);
    unsigned ny = floor((x[1] - y_min)/cell_y_dim);
    unsigned nz = floor((x[2] - z_min)/cell_z_dim);

    /* Check for edge cases */
    if(nx == Nx) { nx = Nx - 1; }
    if(ny == Ny) { ny = Ny - 1; }
    if(nz == Nz) { nz = Nz - 1; }

    /* We know know which bucket this particle goes into. */
    Bucket_Indicies_Body_B[i] = nx + ny*Nx + nz*Nx*Ny;

    /* Update that bucket's B counter */
    #pragma omp atomic
    Buckets[nx + ny*Nx + nz*Nx*Ny].Counter_B++;
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  #ifdef OPERATION_COUNT
    // 3 subractions, 3 divisions per cycle of both loops above.
    #pragma omp single nowait
    {
      #pragma omp atomic update
      OP_Count::Subtraction += 3*(Num_Particles_A + Num_Particles_B);
      #pragma omp atomic update
      OP_Count::Division += 3*(Num_Particles_A + Num_Particles_B);
    } // #pragma omp single
  #endif



  //////////////////////////////////////////////////////////////////////////////
  /* Sort the particles into their corresponding buckets. */
  #pragma omp sections
  {
    #pragma omp section
    {
      /* First, lets set up the Array_A arrays in each bucket (the ith one of
      these holds the indicies of the particles of Body A that are in the ith
      bucket). */
      for(unsigned i = 0; i < Nx*Ny*Nz; i++) { Buckets[i].Array_A.Set_Length(Buckets[i].Counter_A); }

      /* Now, cycle through the elements of Bucket_Indicies_Body_A. Consider
      the ith iteration of this cycle. The ith element of this array tells us
      which bucket particle i of body A belongs in. We, therefore, add its
      index to the corresponding bucket's Array_A array. Importantly, we use the
      Counter_A member of the corresponding bucket to keep track of where we
      should put the particle index in the corresponding bucker's Array_A array
      (remember 0 indexing!) */
      for(unsigned i = 0; i < Num_Particles_A; i++) {
        const unsigned Bucket_Index = Bucket_Indicies_Body_A[i];
        Buckets[Bucket_Index].Counter_A--;
        Buckets[Bucket_Index].Array_A[Buckets[Bucket_Index].Counter_A] = i;
      } // for(unsigned i = 0; i < Num_Particles_A; i++) {
    } // #pragma omp section

    #pragma omp section
    {
      /* First, lets set up the Array_B arrays in each bucket (the ith one of
      these holds the indicies of the particles of body B that are in the ith
      bucket). */
      for(unsigned i = 0; i < Nx*Ny*Nz; i++) { Buckets[i].Array_B.Set_Length(Buckets[i].Counter_B); }

      /* Now, cycle through the elements of Bucket_Indicies_Body_B. Consider
      the ith iteration of this cycle. The ith element of this array tells us
      which bucket particle i of body B belongs in. We, therefore, add its
      index to the corresponding bucket's Array_B array. Importantly, we use the
      Counter_B member of the corresponding bucket to keep track of where we
      should put the particle index in the corresponding bucker's Array_B array
      (remember 0 indexing!) */
      for(unsigned i = 0; i < Num_Particles_B; i++) {
        const unsigned Bucket_Index = Bucket_Indicies_Body_B[i];
        Buckets[Bucket_Index].Counter_B--;
        Buckets[Bucket_Index].Array_B[Buckets[Bucket_Index].Counter_B] = i;
      } // for(unsigned i = 0; i < Num_Particles_B; i++) {
    } // #pragma omp section
  } // #pragma omp sections

  #ifdef CONTACT_MONITOR
    /* Check that each bucket got the right number of particles. This occurs
    precisely when a bucket's Counter's are 0. If they are non-zero, then
    the bucket got too many or too few particles. Since this is stored in the
    global Buckets array, only one thread needs to check this. */
    #pragma omp single
    {
      for(unsigned i = 0; i < Nx*Ny*Nz; i++) {
        if(Buckets[i].Counter_A != 0) { printf("Bucket[%d].Counter_A = %d\n",i, Buckets[i].Counter_A); }
        if(Buckets[i].Counter_B != 0) { printf("Bucket[%d].Counter_B = %d\n",i, Buckets[i].Counter_B); }
      } // for(unsigned i = 0; i < Nx*Ny*Nz; i++) {
    } // #pragma omp single
  #endif



  //////////////////////////////////////////////////////////////////////////////
  /* Finally, let's apply the contact algorithm! */

  /* This function implements Particle-Particle 'contact'.
  If a particle from body A is within Simulation::Contact_Distance of a particle
  in body B, then we apply a contact and friction force between those particles.
  The applied force is in the direction of the line between the two particles'
  centers. */
  const double h = Simulation::Contact_Distance;                               //        : mm
  const double h_squared = h*h;                                                //        : mm^2
  const double Shape_Function_Amp = Body_A.Get_Shape_Function_Amplitude();
  Vector r_ij;                                                                 //        : mm Vector
  Vector Grad_W;                                                               //        : 1/mm^4 Vector
  Vector x_i;                                                                  //        : mm Vector
  Vector F_Contact, F_Friction;                                                //        : N Vector
  Vector Relative_Velocity;                                                    //        : mm/s Vector

  /* Variables to help us apply contact and frictional forces in parallel. */
  bool Contact_Flag = false;
  Vector * Body_B_F_Contact_Local = new Vector[Num_Particles_B];               //        : N Vector
  Vector * Body_B_F_Friction_Local = new Vector[Num_Particles_B];              //        : N Vector

  // First, initialize the local contact/frictional force vectors
  for(unsigned i = 0; i < Num_Particles_B; i++) {
    Body_B_F_Contact_Local[i] = {0,0,0};                                       //        : N Vector
    Body_B_F_Friction_Local[i] = {0,0,0};                                      //        : N Vector
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  // Now loop over the buckets
  #pragma omp for schedule(dynamic) collapse(3)
  for(unsigned kb = 0; kb < Nz; kb++) {
    for(unsigned jb = 0; jb < Ny; jb++) {
      for(unsigned ib = 0; ib < Nx; ib++) {
        /* First, check if body A has any particles in this bucket. If not, then
        move onto the next one. */
        Array<unsigned> & Array_A = Buckets[ib + jb*Nx + kb*Nx*Ny].Array_A;
        const unsigned Num_Particles_A_Bucket = Array_A.Get_Length();
        if(Num_Particles_A_Bucket == 0) { continue; }

        /* In this iteration we will apply the contact algorithm for
        particles of body A in bucket ib + jb*Nx + kb*Nx*Ny. Because of the way
        that we set up the buckets, these particles can only come into contact
        with particles in body B that are in this bucket or a bucket that is
        adjacent to this bucket. Thus, we must cycle through these buckets. We
        need to be careful, however, since there may not be buckets in a
        particular direction (depending on the values of kb, jb, ib) */
        unsigned ib_min = (ib == 0)      ? 0  : ib - 1;
        unsigned jb_min = (jb == 0)      ? 0  : jb - 1;
        unsigned kb_min = (kb == 0)      ? 0  : kb - 1;

        unsigned ib_max = (ib == Nx - 1) ? ib : ib + 1;
        unsigned jb_max = (jb == Ny - 1) ? jb : jb + 1;
        unsigned kb_max = (kb == Nz - 1) ? kb : kb + 1;

        /* Cycle through the buckets that contact can occur in. */
        for(unsigned rb = kb_min; rb <= kb_max; rb++) {
          for(unsigned qb = jb_min; qb <= jb_max; qb++) {
            for(unsigned pb = ib_min; pb <= ib_max; pb++) {
              /* Check if body B has any particles in this bucket. If not, then
              continue. */
              Array<unsigned> & Array_B = Buckets[pb + qb*Nx + rb*Nx*Ny].Array_B;
              const unsigned Num_Particles_B_Bucket = Array_B.Get_Length();
              if(Num_Particles_B_Bucket == 0) { continue; }

              /* Cycle through Body A's particles in bucket ib + jb*Nx + kb*Nx*Ny
              and Body B's particles in bucket ib + jb*Nx + kb*Nx*Ny */
              for(unsigned A_particle_index = 0; A_particle_index < Num_Particles_A_Bucket; A_particle_index++) {
                const unsigned i = Array_A[A_particle_index];

                // Skip broken particles
                if(Body_A[i].Get_D() >= 1) { continue; }

                double V_i = Body_A[i].Get_Volume();                           //        : mm^3
                const double KV_i = K*V_i;                                     //        : N*mm
                x_i = Body_A[i].Get_x();                                       //        : mm Vector

                for(unsigned B_Particle_index = 0; B_Particle_index < Num_Particles_B_Bucket; B_Particle_index++) {
                  const unsigned j = Array_B[B_Particle_index];
                  r_ij = x_i - Body_B[j].Get_x();                              //        : mm Vector

                  // Check if |rij| < h. Note that this is equivalent to rij dot rij < h^2
                  // If so, add contact force
                  if(Dot_Product(r_ij, r_ij) < h_squared) {
                    /* First, set the 'contact flag' to true. This flag is
                    designed to improve perfomance. The idea here is that there
                    will be no contact for most of the time steps. Because of
                    this, it doesn't make sense to run throug the critical for
                    loop (to add each thread's contact and friction forces to
                    the particles). Thus, we only run that critical loop if this
                    flag is true. */
                    Contact_Flag = true;

                    // Calculate the contact force
                    double V_j = Body_B[j].Get_Volume();                       //        : mm^3
                    double Mag_r_ij = Magnitude(r_ij);                         //        : mm
                    double h_minus_Mag_r_ij = h - Mag_r_ij;                    //        : mm
                    Grad_W = (-3.*(Shape_Function_Amp)*(h_minus_Mag_r_ij*h_minus_Mag_r_ij)/Mag_r_ij)*(r_ij);   // 1/mm^4 Vector

                    #ifdef OPERATION_COUNT
                      /* 3 multiplications, 1 division to calculate Grad_W (the
                      final multiplication uses operator overloading and is
                      counted elsewhere. */
                      #pragma omp atomic update
                      OP_Count::Multiplication += 3;
                      #pragma omp atomic update
                      OP_Count::Division += 1;
                    #endif

                    /* Now apply the force to the two interacting bodies (Note
                    the forces are equal and opposite). We have to apply the
                    contact force on B using a critical region. The reson for
                    this is that multiple threads may try to update the same
                    particle's contact force at the same time. Because different
                    thereads operate on different buckets of particles of A, and
                    each particle is in just one bucket, this issue can not
                    happen for particles in body A. However, since different
                    threads can work with the same buckets of B particles, it is
                    possible for two threads to try and update a particle's
                    contact force at the same time, thereby creating a race
                    condition. To fix this, each thread gets its own 'Local
                    contact force array'. each thread stores the force that it
                    would apply to the particles in body B in this array.
                    Once we have cycled through all of body A's particles, we
                    then 'reduce' these force arrays. One by one, the threads
                    add their contributions to body_b's contact forces. This
                    runs about 2x faster than if we had inserted the critical
                    region into this loop. */
                    F_Contact = (KV_i*V_j)*Grad_W;                             //        : N Vector
                    Body_A[i].Force_Contact -= F_Contact;                      //        : N Vector
                    Body_B_F_Contact_Local[j] += F_Contact;                    //        : N Vector

                    /* Now let's calculate the frictional force. To do this,
                    we first need to get the unit vector of the relative
                    velocity between bodies A and B as well as the magnitude of
                    the contact force */
                    double Mag_F_Contact = F_Contact.Magnitude();              //        : N
                    Relative_Velocity = Body_A[i].V - Body_B[j].V;             //        : mm/s Vector

                    // Now we can calculate the frictional force and apply it to the two bodies
                    F_Friction = ((-1*Simulation::Friction_Coefficient*Mag_F_Contact) / Relative_Velocity.Magnitude())*(Relative_Velocity);
                    Body_A[i].Force_Friction += F_Friction;                    //        : N Vector
                    Body_B_F_Friction_Local[j] -= F_Friction;                  //        : N Vector

                    #ifdef OPERATION_COUNT
                      /* F_Contact: 1 multiplication (multiplication by Grad_W uses operator overloading)
                      F_Friction:   2 multiplications, 1 division (multiplycatiom by Relative_Velocity uses operator overloading) */
                      #pragma omp atomic update
                      OP_Count::Multiplication += 3;
                      #pragma omp atomic update
                      OP_Count::Division += 1;
                    #endif
                  } // if(Magnitude(r_ij) < h) {
                } // for(unsigned B_Particle_index = 0; B_Particle_index < Num_Particles_B_Bucket; B_Particle_index++) {
              } // for(unsigned A_particle_index = 0; A_particle_index < Num_Particles_A_Bucket; A_particle_index++) {
            } // for(unsigned pb = ib_min; pb <= ib_max; pb++) {
          } // for(unsigned qb = jb_min; qb <= jb_max; qb++) {
        } // for(unsigned rb = kb_min; rb <= kb_max; rb++) {
      } // for(unsigned ib = 0; ib < Nx; ib++) {
    } // for(unsigned jb = 0; jb < Ny; jb++) {
  } // for(unsigned kb = 0; kb < Nz; kb++) {

  /* Now that we have finished the loop above, each thread has its
  contribution to the contact force for each of Body_B's particles. We therefore
  add these contributions to Body_B's particles one by one (using a critical
  region) */
  #pragma omp critical
  if(Contact_Flag == true) {
    for(unsigned i = 0; i < Num_Particles_B; i++) {
      Body_B[i].Force_Contact += Body_B_F_Contact_Local[i];                    //        : N Vector
      Body_B[i].Force_Friction += Body_B_F_Friction_Local[i];                  //        : N Vector
    } // for(unsigned i = 0; i < Num_Particles_B; i++) {
  } // if(Contact_Flag == true) {

  // Now free any dynamically allocated memory.
  delete [] Body_B_F_Contact_Local;
  delete [] Body_B_F_Friction_Local;

  /* We put a barrier here so that each particle's contact and friction forces
  have been set before deallocating the global arrays and returning. The next
  function in the simulaiton (update_x) will only work correctly if each
  particle's contact and friction forces have been set. */
  #pragma omp barrier

  /* Delete the shared arrays. Why don't we need a barrier? Because there's a
  barrier just before this. The threads can't pass that barrier until they are
  done adding their part of the Contact/Friction forces into the bodies. Thus,
  once they pass this barrier, they are ready to work on update_x. All this
  single directive does is free the dynamic arrays. Thus, the other threads can
  get started on Update_x before this completes */
  #pragma omp single nowait
  {
    delete [] Buckets;
    delete [] Bucket_Indicies_Body_A;
    delete [] Bucket_Indicies_Body_B;
  } // #pragma omp single
} // void Body::Contact_New(Body & Body_A, Body & Body_B) {
