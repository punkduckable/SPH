#include "Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Array.h"
#include "Errors.h"
#include <assert.h>
#include <math.h>

struct Neighbor_Particle_Bucket {
  Array<unsigned> Indices_Array{}; //holds the indices of the particles which are in this bucket
  unsigned Counter = 0;
}; // typdef struct Neighbor_Particle_Bucket

/* Global (shared) variables for the buckets. */
static Neighbor_Particle_Bucket * Buckets;        // Array of buckets
static unsigned * Bucket_Indicies;        // Which bucket each particle of the body goes to.
static double Buffer[6];



////////////////////////////////////////////////////////////////////////////////
// Neighbor methods!

void Body::Set_Neighbors(const unsigned i, const Array<unsigned> & Neighbor_IDs_In) {
  /* This function is used to set the neighbors for each particle in a body. */

  /* First check if this particle already has neighbors. This function should
  only be called if the neighbors have not been set. The reason for this is
  that this method allocates pointers. If the pointers have already been set,
  then allocating them again will cause a memory leak. */
  assert(Particles[i].Neighbors_Are_Set == false);

  // Set Num_Neighbors using input
  unsigned Num_Neighbors_In = Neighbor_IDs_In.Get_Length();
  Particles[i].Num_Neighbors = Num_Neighbors_In;

  /* Next, check that Num_Neighbors_In > 0. if Num_Neighbors_In = 0, then there are no neighbors. */
  if(Num_Neighbors_In == 0) {
    printf("You didn't supply any neighbors for particle %u! Particle %u is being damaged.\n", i, i);
    Particles[i].D = 1;
    return;
  } // if(Num_Neighbors_In == 0) {



  //////////////////////////////////////////////////////////////////////////////
  /* Now that we know our neighbors IDs, we can set Particle[i]'s neighbor list.
  After that has been set, we can call Set_Neighbor_Dependent_Members to set up
  the members of particles[i] that can not be set until the neighbors are
  known (such as r, R, W, and Grad_W, etc...) */

  unsigned * Neighbor_IDs = new unsigned[Num_Neighbors_In];
  for(unsigned j = 0; j < Num_Neighbors_In; j++) { Neighbor_IDs[j] = Neighbor_IDs_In[j]; }
  Particles[i].Neighbor_IDs = Neighbor_IDs;

  Set_Neighbor_Dependent_Members(i);
} // void Body::Set_Neighbors(const unsigned i, const Array<unsigned> & Neighbor_IDs_In) {



void Body::Set_Neighbor_Dependent_Members(const unsigned i) {
  /* Function description:
  This function is designed to set up the members of a particle that can not be
  defined until the particle's members are known (such as R, Mag_R, W, A_Inv,
  etc...).

  This function should NOT be called until the ith particle's neighbor list
  has been set! */
  unsigned Num_Neighbors = Particles[i].Get_Num_Neighbors();

  // Allocate memory for the Dynamic arrays
  Vector * R = new Vector[Num_Neighbors];                                      //        : mm Vector
  double * Mag_R = new double[Num_Neighbors];                                  //        : mm
  double * W = new double[Num_Neighbors];                                      //        : unitless
  Vector * Grad_W = new Vector[Num_Neighbors];                                 //        : 1/mm Vector

  // Declare some variables
  const double h = (*this).Support_Radius;
  int Neighbor_ID;                               // Keep track of current particle
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,                                // Shape Tensor (zero initialized)      : unitless Tensor
           0,0,0,
           0,0,0};

  // Loop through each neighbor, determine relevant information
  for(unsigned j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = Particles[i].Get_Neighbor_IDs(j);        // Get Neighbor ID of the jth neighbor of particle i

    // Calculate displacement vectors
    R[j] = Particles[Neighbor_ID].X - Particles[i].X;      // Reference displacement vector        : mm Vector
    Mag_R[j] = R[j].Magnitude();                           // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amplitude*(h - Mag_R[j])                             //        : unitless
                                   *(h - Mag_R[j])
                                   *(h - Mag_R[j]);

    Grad_W[j] = -3*Shape_Function_Amplitude                                    //        : 1/mm Vector
                  *((h - Mag_R[j])*(h - Mag_R[j]))
                  *(R[j] / Mag_R[j]);

    // Add in the Current Neighbor's contribution to the Shape tensor
    V_j = Particles[Neighbor_ID].Volume;                   // Neighbor Volume            : mm^3
    A += Dyadic_Product((V_j*Grad_W[j]), R[j]);                                //        : unitless Tensor
  } // for(unsigned j = 0; j < Num_Neighbors_In; j++) {


  //////////////////////////////////////////////////////////////////////////////
  // Now set the ith particle's members
  Particles[i].R = R;                                                          //        : mm Vector
  Particles[i].Mag_R = Mag_R;                                                  //        : mm
  Particles[i].W = W;                                                          //        : unitless
  Particles[i].Grad_W = Grad_W;                                                //        : 1/mm Vector
  Particles[i].A_Inv = A^(-1);                                                 //        : unitless Tensor

  // Now that neighbors have been set, we set 'Neighbors_Are_Set' to true
  Particles[i].Neighbors_Are_Set = true;
} // void Body::Set_Neighbor_Dependent_Members(const unsigned i) {



void Body::Find_Neighbors_New(void){
  /* This function finds the neighbors for each particle and stores them in a
  linked list.

  First, we partition the spatial domain. To do this, we first determine the
  maximum and minimum x, y, and z coordinates of the live (damage < 1) particles
  This gives us a box in which the particles of the two bodies live. We
  then divide the x, y, and z dimensions of this box into smaller boxes,
  called cells, each one of which has a side length that is barely greater than
  the support radius. We then allocate an array of buckets with one bucket per
  cell. We then cycle through the particles of A and B, determining which bucket
  each particle belongs to.

  This function assumes that each particle in (*this) has had it's position
  and volume set (these are used to determine neighbors and set neighbor
  dependent members). It also assumes that the body's support radius has
  been set (which is used to find neighbors) */

  //////////////////////////////////////////////////////////////////////////////
  /* Determine how many buckets we need */

  // Get parameters from the body
  const unsigned Num_Particles = (*this).Get_Num_Particles();

  // Declare some variables.
  double x_max, x_min;
  double y_max, y_min;
  double z_max, z_min;

  // Initialize the variables.
  Vector x = (*this)[0].Get_x();
  x_max = x[0];
  x_min = x_max;
  y_max = x[1];
  y_min = y_max;
  z_max = x[2];
  z_min = z_max;

  /* Determine the maximum and minimum x, y, z coordinate of the particles
  in the body. Each thread only searches through a subset of the
  particles. Once this is done, we will "reduce" them together in the
  buffer to determine the global minimum. */
  #pragma omp for nowait
  for(unsigned i = 0; i < Num_Particles; i++) {
    x = (*this)[i].Get_x();

    /* Note: if x[0] > x_max, then we can't also have x[0] < x_min (this relies
    on the fact that x_min <= x_max) */
    if     (x[0] > x_max) { x_max = x[0]; }
    else if(x[0] < x_min) { x_min = x[0]; }

    if     (x[1] > y_max) { y_max = x[1]; }
    else if(x[1] < y_min) { y_min = x[1]; }

    if     (x[2] > z_max) { z_max = x[2]; }
    else if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

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

  #ifdef NEIGHBOR_MONITOR
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
  cell to be >= the support radius. By doing this, particles
  can only be neighbors with particles in their bucket or in buckets
  that are adjacent to their bucket. In general, we want the buckets to be as
  small as possible (so that there are as few particles to check for neighhbors as
  possible). Let's focus on the x coordinate. Let Nx denote the number of
  cell in the x direction. We want Nx to be the largest natural number
  such that Support_Radius <= (x_max - x_min)/Nx. A little thought reveals
  that this occurs precisely when Nx = floor((x_max - x_min)/Support_Radius).
  A similar result holds for the y and z directions. */
  const double h = (*this).Support_Radius; //set h to the support radius for convenience
  const unsigned Nx = floor((x_max - x_min)/h);
  const unsigned Ny = floor((y_max - y_min)/h);
  const unsigned Nz = floor((z_max - z_min)/h);

  /* Now that we know the number of buckets in each direction, one thread can
  allocate the global arrays */
  #pragma omp single
  {
    Buckets = new Neighbor_Particle_Bucket[Nx*Ny*Nz];
    Bucket_Indicies = new unsigned[Num_Particles];
  } // #pragma omp single

  #ifdef NEIGHBOR_MONITOR
    /* Check that Nx, Ny, Nz were set up properly and that all threads see the
    same dynamically allocated global arrays */
    #pragma omp critical
    {
      printf("(Nx, Ny, Nz) = (%d, %d, %d)\n", Nx, Ny, Nz);
      printf("Buckets = %p\n", Buckets);
      printf("Bucket_Indicies = %p\n", Bucket_Indicies);
    } // #pragma omp critical
  #endif

  //////////////////////////////////////////////////////////////////////////////
  /* Now, determine which bucket each particle belongs to. */

  /* First, we calculate some variables (see the next comment for an
  explanation) */
  const double cell_x_dim = (x_max - x_min)/Nx;
  const double cell_y_dim = (y_max - y_min)/Ny;
  const double cell_z_dim = (z_max - z_min)/Nz;

  #pragma omp for
  for(unsigned i = 0; i < Num_Particles; i++) {
    /* First, we need to determine which bucket our particle belongs in. Let's
    focus on the x coordinate. Each bucket has an x-dimension length of
    (x_max - x_min)/Nx, which we call cell_x_dim. The x coordinates of
    the bucket for the ith particle is the number of units of length
    cell_x_dim that fit between x_min and the x coordinate of the particle
    (think about it). A little thought reveals that this is precisely
    floor((Particle[i].x[0] - x_min)/cell_x_dim)
    A similar result holds for y and z. */
    x = (*this).Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/cell_x_dim);
    unsigned ny = floor((x[1] - y_min)/cell_y_dim);
    unsigned nz = floor((x[2] - z_min)/cell_z_dim);

    /* We run into a bit of a problem if a particle's x coordinate is equal to
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
    Bucket_Indicies[i] = nx + ny*Nx + nz*Nx*Ny;

    /* Update that bucket's counter */
    #pragma omp atomic
      Buckets[nx + ny*Nx + nz*Nx*Ny].Counter++;
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  //////////////////////////////////////////////////////////////////////////////
  /* Sort the particles into their corresponding buckets. */
  #pragma omp single
  {
    /* First, lets set up the Indices_Array arrays in each bucket (the ith one of
    these holds the indicies of the particles of Body A that are in the ith
    bucket). */
    for(unsigned i = 0; i < Nx*Ny*Nz; i++) { Buckets[i].Indices_Array.Set_Length(Buckets[i].Counter); }

    /* Now, cycle through the elements of Bucket_Indicies. Consider
    the ith iteration of this cycle. The ith element of this array tells us
    which bucket particle i of body A belongs in. We, therefore, add its
    index to the corresponding bucket's Indices_Array array. Importantly, we use the
    Counter member of the corresponding bucket to keep track of where we
    should put the particle index in the corresponding bucket's Indices_Array array
    (remember 0 indexing!) */
    for(unsigned i = 0; i < Num_Particles; i++) {
      const unsigned Bucket_Index = Bucket_Indicies[i];
      Buckets[Bucket_Index].Counter--;
      Buckets[Bucket_Index].Indices_Array[Buckets[Bucket_Index].Counter] = i;
    } // for(unsigned i = 0; i < Num_Particles; i++) {
  } // #pragma omp single

  #ifdef NEIGHBOR_MONITOR
    /* Check that each bucket got the right number of particles. This occurs
    precisely when a bucket's Counter's are 0. If they are non-zero, then
    the bucket got too many or too few particles. Since this is stored in the
    global Buckets array, only one thread needs to check this. */
    #pragma omp single
    {
      for(unsigned i = 0; i < Nx*Ny*Nz; i++) {
        if(Buckets[i].Counter != 0) { printf("Bucket[%d].Counter = %d\n",i, Buckets[i].Counter); }
      } // for(unsigned i = 0; i < Nx*Ny*Nz; i++) {
    } // #pragma omp single
  #endif

  //////////////////////////////////////////////////////////////////////////////
  /* Finally, let's find the neighbors! */

  /* A particle is considered to be neighbors with another particle if they are
  in the same body and the distance between them is less than the support radius.*/
  const double h_squared = h*h;
  Vector Rij;
  List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors

  //loop over the buckets
  #pragma omp for schedule(dynamic) collapse(3)
  for(unsigned kb = 0; kb < Nz; kb++) {
    for(unsigned jb = 0; jb < Ny; jb++) {
      for(unsigned ib = 0; ib < Nx; ib++) {
        /* First, check if body A has any particles in this bucket. If not, then
        move onto the next one. */

        Array<unsigned> & Indices_Array = Buckets[ib + jb*Nx + kb*Nx*Ny].Indices_Array;
        const unsigned Num_Particles = Indices_Array.Get_Length();
        if(Num_Particles == 0) { continue; }

        /* In this iteration we will find the neighbors for
        particles in bucket ib + jb*Nx + kb*Nx*Ny. Because of the way
        that we set up the buckets, these particles can only be neighbors with
        particles that are in this bucket or a bucket that is
        adjacent to this bucket. Thus, we must cycle through these buckets. We
        need to be careful, however, since there may not be buckets in a
        particular direction (depending on the values of kb, jb, ib) */
        unsigned ib_min = (ib == 0)      ? 0  : ib - 1;
        unsigned jb_min = (jb == 0)      ? 0  : jb - 1;
        unsigned kb_min = (kb == 0)      ? 0  : kb - 1;

        unsigned ib_max = (ib == Nx - 1) ? ib : ib + 1;
        unsigned jb_max = (jb == Ny - 1) ? jb : jb + 1;
        unsigned kb_max = (kb == Nz - 1) ? kb : kb + 1;

        /* Cycle through particles in main bucket */
        for(unsigned int Main_particle_index = 0; Main_particle_index < Num_Particles; Main_particle_index++) {
          unsigned int i = Indices_Array[Main_particle_index]; //index of particle whose neighbors we are trying to find

          /* Cycle through the buckets that neighbors may be in. */
          for(unsigned rb = kb_min; rb <= kb_max; rb++) {
            for(unsigned qb = jb_min; qb <= jb_max; qb++) {
              for(unsigned pb = ib_min; pb <= ib_max; pb++) {

                Array<unsigned> & Neighbors_Array = Buckets[pb + qb*Nx + rb*Nx*Ny].Indices_Array;
                const unsigned Num_Particles_Neighbor_Bucket = Neighbors_Array.Get_Length();
                if(Num_Particles_Neighbor_Bucket == 0) { continue; }

                /*cycle through particles in Neighbors_Array (corresponding to bucket pb + qb*Nx + rb*Nx*Ny) */
                for(unsigned int Neighbor_particle_index = 0; Neighbor_particle_index < Num_Particles_Neighbor_Bucket; Neighbor_particle_index++) {
                  unsigned int j = Neighbors_Array[Neighbor_particle_index]; //index of neighbor particle

                  //make sure that i and j don't correspond to the same particle
                  //(a particle can't be its own neighbor)
                  if (i != j) {
                    Rij = (*this).Particles[i].Get_X() - (*this).Particles[j].Get_X();
                    /* If |Rij| < h then i and j are neighbors */
                    if (h_squared > Dot_Product(Rij, Rij)){Particle_Neighbor_List.Push_Back(j); }
                  } // if (i != j) {
                } //for(unsigned int Neighbor_particle_index = 0; Neighbor_particle_index < Num_Particles_Neighbor_Bucket; Neighbor_particle_index++)
              } //for(unsigned pb = ib_min; pb <= ib_max; pb++)
            } //for(unsigned qb = jb_min; qb <= jb_max; qb++)
          } //for(unsigned rb = kb_min; rb <= kb_max; rb++)

          /* Now that we have the neighbor ID list, we can make it into an array.
          This is done using the Array class' list constructor. This will also empty
          the Particle_Neighbor_list so that it is ready for the next i. See Array.h */
          Array<unsigned> Neighbor_IDs(Particle_Neighbor_List);

          // Now send the Neighbor list to the particle
          Set_Neighbors(i, Neighbor_IDs);
        } // for(unsigned int Main_particle_index = 0; Main_particle_index < Num_Particles; Main_particle_index++) {
      } //for(unsigned ib = 0; ib < Nx; ib++)
    } //for(unsigned jb = 0; jb < Ny; jb++)
  } //for(unsigned kb = 0; kb < Nz; kb++)
} // void Body::Find_Neighbors_New(void){



bool Body::Are_Neighbors(const unsigned i, const unsigned j) const {
  /* This function checks if h > |Rj|. Here, Rj is simply the displacement of
  particle i relative to particle j: Rj = Xj - Xi. Xj = P1.X, Xi = P2.X. if
  h > |Rj| then P1 and P2 are in each other's support radius, so P1 is a
  neighbor of P2.

  It should be noted that, since h and |Rj| are positive numbers, if h>|Rj|
  then h^2>|Rj|^2. We can compute this second condition using a dot product
  (which is far easier than finding the magnitude) */

  /* A particle can not be its own neighbor. If i = j, then something went wrong.
  and we should abort. */
  assert(i != j);

  const double h = (*this).Support_Radius;
  const Vector Rj = (*this).Particles[i].Get_X() - (*this).Particles[j].Get_X();
  return ( h*h > Dot_Product(Rj, Rj) );
} // bool Body::Are_Neighbors(const unsigned i, const unsigned j) const {


//old Find_Neighbors function
void Body::Find_Neighbors(void) {
  /* This function finds the neighbors for each particle in the (*this) body.
  It assumes that the particles in (*this) have had their positions set. */

  unsigned i,j;                              // Loop index variables
  List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors

  // Cycle through the particles
  #pragma omp for
  for(i = 0; i < Num_Particles; i++) {

    // For each particle, cycle through the potential neighbors (every particle)
    for(j = 0; j < Num_Particles; j++) {
      // ith particle is not its own neighbor.
      if(j == i) { continue; }

      // Test if jth particle is inside support radius of ith particle. If so,
      // add P_j to P_i's neighbor list.
      if(Are_Neighbors(i, j)) { Particle_Neighbor_List.Push_Back(j); }
    } // for(unsigned j = 0; j < Num_Particles; j++) {

    // Now that we have the neighbor ID list, we can make it into an array.
    // This is done using the Array class' list constructor. See Array.h
    Array<unsigned> Neighbor_IDs(Particle_Neighbor_List);

    // Now sent the Neighbor list to the particle
    Set_Neighbors(i, Neighbor_IDs);
  } // for(unsigned i = 0; i < Num_Particles; i++) {
} // void Body::Find_Neighbors(void) {



void Body::Find_Neighbors_Box(void) {
  /* This function is a modified version of the Neighbor List generating
  function that is specialized for Box particle geometries.

  Let us establish a few definitions:
  By a 'Vertical Layer' we mean a sheet of particles that all have the same x
  coordinate
  By a 'Vertical column' we mean a set of particles with the same x and z
  coordinates.
  By a 'Row' we mean a set of particles with the same x and y coordinates.

  This function assumes that the particles are stored in 'Vertical-Column'
  major, 'Vertical Layer' semi-major order. This menas that vertical columns of
  particles are stored in contiguous memory and that vertical columns in the
  same vertical layer are stored in contiguous memory.

  Thus, if working with a cube with sidelength N, the (1,1,1) particle will be
  N*N particles away from the (2,1,1) particle in the, N particles away from the
  (1,1,2) particle and 1 particle away from the (1,2,1) particle in the body.

  So why does this function exist?
  A generic neighbor search is slow. For a given particle to find its neighbors,
  it has no choice but to search through EVERY other particle in the body.
  If the body has M particles then there are a total of M*M
  neighbor tests performed in all. This is highly inefficient. However, if
  we're working with a Box of particles, then the particles are stored in
  a regular grid pattern. Rather than searching through every particle in
  the body, we can just search through the grid elements that are
  close to the current particle! This reduces the number of searches with a M
  body from M*M to M*(d_max^3), where d_max = floor(Support_Radius/IPS)

  So how do you use this function?
  This function is used just like the Generage_Neighbor_List function. We can
  determine the Box dimensions from the Body. These dimensions are
  stored in X_SIDE_LENGTH, Y_SIDE_LENGTH, and Z_SIDE_LENGTH. If the Box has p
  particles in a vertical column then Z_SIDE_LENGTH is p. For a 100x50x200
  Box of particles, X_SIDE_LENGTH is 100, Y_SIDE_LENGTH is 50, and
  Z_SIDE_LENGTH is 200 */

  #pragma omp single
  {
    if((*this).Is_Box == false) {
      char Buf[500];
      sprintf(Buf,
              "Not A Box Exception: thrown by Body::Find_Neighbors_Box\n"
              "Body %s tried to use this function, but %s is not a box! This function\n"
              "can only be called by boxes!\n",
              (*this).Name.c_str(), (*this).Name.c_str());
      throw Not_A_Box(Buf);
    } // if((*this).Is_Box == false) {
  } // #pragma omp single

  const unsigned d_max = floor((*this).Support_Radius / (*this).Inter_Particle_Spacing);

  #pragma omp for collapse(3)
  for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {
    for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
      for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
        unsigned p,q,r;                             // Loop index variables
        unsigned p_min, p_max, q_min, q_max, r_min, r_max;
        List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors

        /* If we are near the upper or lower bound for any one of the 3
        coordinates, then we need to adjust the upper or lower bound of our
        search.

        To understand why we need to do this, suppose that our i coordinate is 3.
        Then, the only smaller i coordinates are 0, 1, and 2. In general, if our
        i coordinate is n then there are n smaller i coordinate values. If the
        support radius is < n, then we need to check the n-d_max i
        coordinates with a smaller i coordinate. Otherwise, we need to check
        all n. */

        // i index (x coordinate) checks
        if(i < d_max) { p_min = 0; }
        else{ p_min = i -  d_max; }

        if(i + d_max > (X_SIDE_LENGTH - 1)) { p_max = X_SIDE_LENGTH - 1; }
        else { p_max = i + d_max; }

        // j index (y coordinate) checks
        if(j < d_max) { q_min = 0; }
        else { q_min = j - d_max; }

        if(j + d_max > (Y_SIDE_LENGTH - 1)) { q_max = Y_SIDE_LENGTH - 1; }
        else { q_max = j + d_max; }

        // k index (z coordinate) checks
        if(k < d_max) { r_min = 0; }
        else { r_min = k - d_max; }

        if(k + d_max > (Z_SIDE_LENGTH - 1)) { r_max = Z_SIDE_LENGTH - 1; }
        else { r_max = k + d_max; }

        // Loop through potential neighbors, generate neighbor list
        for(p = p_min; p <= p_max; p++) {
          for(r = r_min; r <= r_max; r++) {
            for(q = q_min; q <= q_max; q++) {
              // a given particle is NOT its own neighbor
              if(i == p && j == q && k == r) { continue; }

              if(Are_Neighbors(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*(Y_SIDE_LENGTH) + j,
                               p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q))
                Particle_Neighbor_List.Push_Back(p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q);
            } // for(q = q_min; q <= q_max; q++) {
          } // for(r = r_min; r <= r_max; r++) {
        } // for(p = p_min; p <= p_max; p++) {

        /* Now that we have the neighbor list, we can make it into an array. To do
        this, I use the Array class' list constructor. See Array.h */
        Array<unsigned> Neighbor_IDs(Particle_Neighbor_List);

        // Now sent the Neighbor list to the particle
        Set_Neighbors(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*(Y_SIDE_LENGTH) + j, Neighbor_IDs);
      } // for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {
} // void Body::Find_Neighbors_Box(void) {



void Body::Remove_Neighbor(const unsigned i, const unsigned Remove_Neighbor_ID) {
  /* This function is used to remove a neighbor from a particle. */

  // To be able to remove a neighbor, the particle in question must have neighbors!
  if(Particles[i].Neighbors_Are_Set == false || Particles[i].Num_Neighbors == 0) {
    printf("Particle %d has no neighbors! We can't remove %d\n", i, Remove_Neighbor_ID);
    return;
  } // if(P_In.Neighbors_Are_Set == false || P_In.Num_Neighbors == 0) {

  /* Note: We use the term 'Neighbor Arrays' to refer to the dynamic particle
  member varialbes that hold neighbor information (R, W, Grad_W, etc...)

  So how do we remove a neighbor? Simple, we keep every old neighbor except
  for the specified one. We do this by allocating new arrays with one fewer than
  the old number of neighbors! We then cycle through this particles old
  neihbors. For each old neighbor, we check if its ID matches Remove_Neighbor_ID.
  If it doesn't match, then we copy that neighbor's data from the old Neighbor
  arrays into the new Neighbor arrays. If it does match, then we skip that
  neighbor (don't copy its info over). Once we are finished, we delete the old
  Neighbor arrays and point this body's to the new Neighbor arrays. */

  unsigned Num_Neighbors = Particles[i].Num_Neighbors;
  unsigned p = -1;                               // new neighbors index variables
  double Vol_p;                                  // Volume of pth neighbor

  // New neighbor arrays
  unsigned *New_Neighbor_IDs = new unsigned[Num_Neighbors - 1];                //        : unitless
  Vector *New_R = new Vector[Num_Neighbors - 1];                               //        : mm Vector
  double *New_Mag_R = new double[Num_Neighbors - 1];                           //        : mm
  double *New_W = new double[Num_Neighbors - 1];                               //        : 1/(mm^3)
  Vector *New_Grad_W = new Vector[Num_Neighbors - 1];                          //        : 1/(mm^4) Vector
  Tensor New_A{0,0,0,
               0,0,0,
               0,0,0};

  for(unsigned j = 0; j < Num_Neighbors; j++) {
    // Check if jth neighbor ID matches Remove_Neighbor_ID
    if(Particles[i].Neighbor_IDs[j] == Remove_Neighbor_ID) { continue; }

    // If not, then this is a new neighbor, increment j.
    p++;

    // Check if p == Num_Neighbors - 1. If this is the case, then Remove_Particle_ID
    // was NOT one of this particle's neighbors!
    if(p == Num_Neighbors - 1) {
      printf("particle %d was not a neighbor of particle %d\n",Remove_Neighbor_ID, i);
      return;
    } // if(p == Num_Neighbors - 1) {

    // Copy old Neighbor data to New Neighbor arrays.
    New_Neighbor_IDs[p] = Particles[i].Neighbor_IDs[j];                        //        : unitless
    New_R[p]            = Particles[i].R[j];                                   //        : mm Vector
    New_Mag_R[p]        = Particles[i].Mag_R[j];                               //        : mm
    New_W[p]            = Particles[i].W[j];                                   //        : 1/(mm^3)
    New_Grad_W[p]       = Particles[i].Grad_W[j];                              //        : 1/(mm^4) Vector

    // Calculate New shape tensor.
    Vol_p = Particles[Particles[i].Neighbor_IDs[p]].Volume;                    //        : mm^3
    New_A += Dyadic_Product((Vol_p*New_Grad_W[p]), New_R[p]);                  // New shape tensor : unitless Tensor
  } // for(unsigned j = 0; j < Num_Neighbors; j++) {

  // Now that we have our new neighbor arrays, we can replace/delete the old
  // neighbor arrays
  delete [] Particles[i].Neighbor_IDs;
  delete [] Particles[i].R;
  delete [] Particles[i].Mag_R;
  delete [] Particles[i].W;
  delete [] Particles[i].Grad_W;

  Particles[i].Neighbor_IDs = New_Neighbor_IDs;                                //        : unitless
  Particles[i].R = New_R;                                                      //        : mm Vector
  Particles[i].Mag_R = New_Mag_R;                                              //        : mm
  Particles[i].W = New_W;                                                      //        : 1/(mm^3)
  Particles[i].Grad_W = New_Grad_W;                                            //        : 1/(mm^4) Vector

  // Decrement number of neighbors by 1.
  Particles[i].Num_Neighbors--;

  // Now we can calculate the new A^(-1) from New_A.
  Particles[i].A_Inv = New_A^(-1);                                             //        : unitless Tensor
} // void Body::Remove_Neighbor(const unsigned i, const unsigned Remove_Neighbor_ID) {
