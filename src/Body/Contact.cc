#include "Body.h"
#include "Simulation/Simulation.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Array.h"
#include <math.h>

struct Contact_Particle_Bucket {
  List<unsigned> List_A{};
  List<unsigned> List_B{};

  Array<unsigned> Array_A{};
  Array<unsigned> Array_B{};
}; // typdef struct Particle_Bucket {


void Body::Contact(Body & Body_A, Body & Body_B) {
  /* This function calculates the contact force between the particles in both A
  and those in body B.

  First, we partition the spatial domain. To do this, we first determine the
  maximum and minimum x, y, and z coordinates of the live (damage < 1) particles
  This gives us a cuboid in which the particles of the two bodies live. We
  then divide the x, y, and z dimensions of this cuboid into smaller cuboids,
  each one of which has a side length that is barely greater than the contact
  distance. We then allocate an array of buckets with one bucket per sub-cuboid.
  We then cycle through the particles of A and B, sorting them into their
  corresponding buckets (specifically the lists in those buckets).

  Once every particle has been sorted, the bucket lists are used to generate
  arrays in each bucket (we use lists at first because we don't know ahead of
  time how many particles will go into each bucket. We convert these lists into
  arrays so that we can access them quickly).

  We then apply the contact algorithm. The way that these buckets are
  set up, a particle can only contact particles that are in its bucket or a
  bucket which is adjacent to its bucket (27 total buckets to check). This
  greatly reduces the number of particles that we need to check against, thereby
  eliminating needless calculations. */



  //////////////////////////////////////////////////////////////////////////////
  /* Determine how many buckets we need */

  // Get paramaters from Body_A, Body_B
  const unsigned Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned Num_Particles_B = Body_B.Get_Num_Particles();
  const double Shape_Function_Amp = Body_A.Get_Shape_Function_Amplitude();

  // Set up variables.
  double x_max, x_min;
  double y_max, y_min;
  double z_max, z_min;

  Vector x = Body_A.Particles[0].Get_x();
  x_max = x[0];
  x_min = x_max;
  y_max = x[1];
  y_min = y_max;
  z_max = x[2];
  z_min = z_max;

  /* Determine the maximum and minimum x, y, z coordinate of the particles
  in the two bodies. */
  for(unsigned i = 0; i < Num_Particles_A; i++) {
    x = Body_A.Particles[i].Get_x();

    if(x[0] > x_max) { x_max = x[0]; }
    if(x[0] < x_min) { x_min = x[0]; }

    if(x[1] > y_max) { y_max = x[1]; }
    if(x[1] < y_min) { y_min = x[1]; }

    if(x[2] > z_max) { z_max = x[2]; }
    if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  for(unsigned i = 0; i < Num_Particles_B; i++) {
    x = Body_B.Particles[i].Get_x();

    if(x[0] > x_max) { x_max = x[0]; }
    if(x[0] < x_min) { x_min = x[0]; }

    if(x[1] > y_max) { y_max = x[1]; }
    if(x[1] < y_min) { y_min = x[1]; }

    if(x[2] > z_max) { z_max = x[2]; }
    if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  /* We want the dimension (in all three coordinate directions) of the
  sub-cuboids to be >= the contact distance. By doing this, particles
  can only compe into contact with particles in their bucket or in buckets
  that are adjacent to their bucket. In general, we want the buckets to be as
  small as possible (so that there are as few particles to check for contact as
  possible). Let's focus on the x coordinate. Let Nx denote the number of
  sub-cuboids in the x direction. We want Nx to be the largest natural number
  such that Contact_Distance <= (x_max - x_min)/Nx. A little though reveals
  that this occurs precisely when Nx = floor((x_max - x_min)/Contact_Distance).
  A similar result holds for the y and z directions. */
  const unsigned Nx = floor((x_max - x_min)/Simulation::Contact_Distance);
  const unsigned Ny = floor((y_max - y_min)/Simulation::Contact_Distance);
  const unsigned Nz = floor((z_max - z_min)/Simulation::Contact_Distance);

  /* Now that we know the number of buckets in each direction, we can allocate
  the buckets array! */
  Array<struct Contact_Particle_Bucket> Buckets{Nx*Ny*Nz};



  //////////////////////////////////////////////////////////////////////////////
  /* Now, sort the particles of A and B into their corresponding buckets. */

  /* First, we calculate some variables (see the next comment for an
  explanation) */
  const double sub_cuboid_x_dim = (x_max - x_min)/Nx;
  const double sub_cuboid_y_dim = (y_max - y_min)/Ny;
  const double sub_cuboid_z_dim = (y_max - z_min)/Nz;

  for(unsigned i = 0; i < Num_Particles_A; i++) {
    /* First, we need to determine which bucket our particle belongs in. Let's
    focus on the x coordinate. Each bucket has an x-dimension length of
    (x_max - x_min)/Nx, which we call sub_cuboid_x_dim. The x coordinates of
    the bucket for the ith particle is the number of units of length
    sub_cuboid_x_dim that fit between x_min and the x coordinate of the particle
    (think about it). A little thought reveals that this is precisely
    floor((Particle[i].x[0] - x_min)/((x_max - x_min)/Nx))
    A similar result holds for y and z. */
    x = Body_A.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/sub_cuboid_x_dim);
    unsigned ny = floor((x[1] - y_min)/sub_cuboid_y_dim);
    unsigned nz = floor((x[2] - z_min)/sub_cuboid_z_dim);

    Buckets[nx + ny*Nx + nz*Nx*Ny].List_A.Push_Front(i);
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  for(unsigned i = 0; i < Num_Particles_B; i++) {
    x = Body_B.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/sub_cuboid_x_dim);
    unsigned ny = floor((x[1] - y_min)/sub_cuboid_y_dim);
    unsigned nz = floor((x[2] - z_min)/sub_cuboid_z_dim);

    Buckets[nx + ny*Nx + nz*Nx*Ny].List_B.Push_Front(i);
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  /* Now, let's convert the buckets lists into arrays. */
  for(unsigned r = 0; r < Nz; r++) {
    for(unsigned q = 0; q < Ny; q++) {
      for(unsigned p = 0; p < Nx; p++) {
        Buckets[p + q*Nx + r*Nx*Ny].Array_A.Setup_From_List(Buckets[p + q*Nx + r*Nx*Ny].List_A);
        Buckets[p + q*Nx + r*Nx*Ny].Array_B.Setup_From_List(Buckets[p + q*Nx + r*Nx*Ny].List_B);
      } // for(unsigned p = 0; p < Nx; p++) {
    } // for(unsigned q = 0; q < Ny; q++) {
  } // for(unsigned r = 0; r < Nz; r++) {



  //////////////////////////////////////////////////////////////////////////////
  /* Finally, let's apply the contact algorithm! */

  /* This function implements Particle-Particle 'contact'.
  If a particle from body A is within Simulation::Contact_Distance of a particle
  in body B, then we apply a contact and friction force between those particeles.
  The applied force is in the direction of the line between the two particles'
  centers. */
  const double h = Simulation::Contact_Distance;                               //        : mm
  const double h_squared = h*h;                                                //        : mm^2
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
              if(Num_Particles_B == 0) { continue; }

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
                    Grad_W = (-3*(Shape_Function_Amp)*(h_minus_Mag_r_ij*h_minus_Mag_r_ij)/Mag_r_ij)*(r_ij);   // 1/mm^4 Vector

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
    } // for(unsignd i = 0; i < Num_Particles_B; i++) {
  } // if(Contact_Flag == true) {

  // Now free any dynamically allocated memory.
  delete [] Body_B_F_Contact_Local;
  delete [] Body_B_F_Friction_Local;

  /* We put a barrier here so that each particle's contact and friction forces
  have been set before returning. The next function in the simulaiton (update_x)
  will only work correctly if each particle's contact and friction forces have
  been set. */
  #pragma omp barrier
} // void Body::Contact_New(Body & Body_A, Body & Body_B) {
