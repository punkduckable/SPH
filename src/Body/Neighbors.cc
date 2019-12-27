#include "Body.h"

////////////////////////////////////////////////////////////////////////////////
// Neighbor methods!

void Body::Set_Neighbors(const unsigned i, const unsigned Num_Neighbors_In, const unsigned * Neighbor_ID_Array) {
  /* First check if this particle already has neighbors. This function should
  only be called if the neighbors have not been set. The reason for this is
  that this method allocates pointers. If the pointers have already been set,
  then allocating them again will cause a memory leak. */
  assert(Particles[i].Neighbors_Are_Set == false);

  // Set Num_Neighbors using input
  Particles[i].Num_Neighbors = Num_Neighbors_In;

  /* Next, check that Num_Neighbors_In > 0. if Num_Neighbors_In = 0, then there are no neighbors. */
  if(Num_Neighbors_In == 0) {
    printf("You didn't supply any neighbors for particle %u! I'm damaging this particle.\n", i);
    Particles[i].D = 1;
    return;
  } // if(Num_Neighbors_In == 0) {



  //////////////////////////////////////////////////////////////////////////////
  /* Now that we know our neighbors IDs, we can figure out everything that we
  want to know about them. We an set the Neighbor_IDs, r, R, W, and Grad_W
  members. These can be used to calculate the shape matrix (and its inverse)! */


  // Allocate memory for the Dynamic arrays
  unsigned * Neighbor_IDs = new unsigned[Num_Neighbors_In];
  Vector * R = new Vector[Num_Neighbors_In];                                   //        : mm Vector
  double * Mag_R = new double[Num_Neighbors_In];                               //        : mm
  double * W = new double[Num_Neighbors_In];                                   //        : unitless
  Vector * Grad_W = new Vector[Num_Neighbors_In];                              //        : 1/mm Vector

  // Allocate some variables
  int Neighbor_ID;                               // Keep track of current particle
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,                                // Shape Tensor (zero initialized)      : unitless Tensor
           0,0,0,
           0,0,0};

  // Loop through each neighbor, determine relevant information
  for(unsigned j = 0; j < Num_Neighbors_In; j++) {
    Neighbor_ID = Neighbor_ID_Array[j];          // Get Neighbor ID (index in Particles array)
    Neighbor_IDs[j] = Neighbor_ID;               // Set jth element of Neighbor_IDs member

    // Calculate displacement vectors
    R[j] = Particles[Neighbor_ID].X - Particles[i].X;      // Reference displacement vector        : mm Vector
    Mag_R[j] = R[j].Magnitude();                 // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amplitude*(h - Mag_R[j])                             //        : unitless
                                   *(h - Mag_R[j])
                                   *(h - Mag_R[j]);
    Grad_W[j] = -3*Shape_Function_Amplitude                                    //        : 1/mm Vector
                  *((h - Mag_R[j])*(h - Mag_R[j]))
                  *(R[j] / Mag_R[j]);

    // Add in the Current Neighbor's contribution to the Shape tensor
    V_j = Particles[Neighbor_ID].Vol;            // Neighbor Volume                      : mm^3
    A += Dyadic_Product((V_j*Grad_W[j]), R[j]);                                //        : unitless Tensor
  } // for(unsigned j = 0; j < Num_Neighbors_In; j++) {



  //////////////////////////////////////////////////////////////////////////////
  // Now set the ith particle's members
  Particles[i].Neighbor_IDs = Neighbor_IDs;
  Particles[i].R = R;                                                          //        : mm Vector
  Particles[i].Mag_R = Mag_R;                                                  //        : mm
  Particles[i].W = W;                                                          //        : unitless
  Particles[i].Grad_W = Grad_W;                                                //        : 1/mm Vector
  Particles[i].A_Inv = A^(-1);                                                 //        : unitless Tensor


  // Now that neighbors have been set, we set 'Neighbors_Are_Set' to true
  Particles[i].Neighbors_Are_Set = true;
} // void Body::Set_Neighbors(const unsigned i, const unsigned Num_Neighbors, const unsigned * Neighbor_ID_Array) {



bool Body::Are_Neighbors(const unsigned i, const unsigned j) const {
  /* This function checks if h > |Rj|. Here, Rj is simply the displacement of
  particle i relative to particle j: Rj = Xj - Xi. Xj = P1.X, Xi = P2.X. if
  h > |Rj| then P1 and P2 are in each other's support radius, so P1 is a
  neighbor of P2.

  It should be noted that, since h and |Rj| are positive numbers, if h>|Rj|
  then h^2>|Rj|^2. We can compute this second condition using a dot product
  (which is far easier than finding the magnitude)*/

  const Vector Rj = (*this).Particles[i].Get_X() - (*this).Particles[j].Get_X();
  return ( h*h > Dot_Product(Rj, Rj) );
} // bool Body::Are_Neighbors(const unsigned i, const unsigned j) const {



void Body::Find_Neighbors(void) {
  unsigned i,j;                              // Loop index variables
  List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors
  unsigned *Neighbor_IDs;                    // Array that holds final list of neighbors

  // Cycle through the particles
  for(i = 0; i < Num_Particles; i++) {

    /* For each particle, cycle through the potential neighbors (every particle) */
    for(j = 0; j < Num_Particles; j++) {
      // ith particle is not its own neighbor.
      if(j == i) { continue; }

      // Test if jth particle is inside support radius of ith particle. If so,
      // add P_j to P_i's neighbor list.
      if(Are_Neighbors(i, j)) { Particle_Neighbor_List.Add_Back(j); }
    } // for(unsigned j = 0; j < Num_Particles; j++) {

    /* Now that we have the neighbor list, we can make it into an array. To do
    this, we allocate an array whose length is equal to the length of the
    neighbor list. We then populate this array with the elements of the list
    and finally send this off to the particle (whose neighbors we found) */
    unsigned Num_Neighbors = Particle_Neighbor_List.Node_Count();
    Neighbor_IDs = new unsigned[Num_Neighbors];

    for(j = 0; j < Num_Neighbors; j++) { Neighbor_IDs[j] = Particle_Neighbor_List.Remove_Front(); }

    // Now sent the Neighbor list to the particle
    Set_Neighbors(i, Num_Neighbors, Neighbor_IDs);

    /* Now free Neighbor_IDs array for next particle! */
    delete [] Neighbor_IDs;
  } // for(unsigned i = 0; i < Num_Particles; i++) {
} // void Body::Find_Neighbors(void) {



void Body::Find_Neighbors_Cuboid(void) {
  assert(Is_Cuboid == true);

  /* This function is a modified version of the Neighbor List generating
  function that is specialized for Cuboid particle geometries.

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
  N*N particles away from the (2,1,1) partilce in the, N particles away from the
  (1,1,2) particle and 1 particle away from the (1,2,1) particle in the
  Particles array.

  So why does this function exist?
  A generic neighbor search is slow. For a given particle to find its neighbors,
  it has no choice but to search through EVERY other particle in the Particles
  array. If the Particles array has M particles then there are a total of M*M
  neighbor tests performed in all. This is highly inefficient. However, if
  we're working with a cuboid of particles, then the particles are stored in
  a regular grid pattern. Rather than searching through every particle in
  the Particles array, we can just search through the grid elements that are
  close to the current particle! This reduces the number of searches with a M
  particle array from M*M to M*(Support_Radius^3), where Support radius is in
  units of inter particle spacings.

  So how do you use this function?
  This function is used just like the Generage_Neighbor_List function. We can
  determine the cuboid dimensions from the Body. These dimensions are
  stored in X_SIDE_LENGTH, Y_SIDE_LENGTH, and Z_SIDE_LENGTH. If the cuboid has p
  particles in a vertical column then Z_SIDE_LENGTH is p. For a 100x50x200
  cuboid of particles, X_SIDE_LENGTH is 100, Y_SIDE_LENGTH is 50, and
  Z_SIDE_LENGTH is 200 */

  for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {
    for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
      for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
        unsigned p,q,r;                             // Loop index variables
        unsigned p_min, p_max, q_min, q_max, r_min, r_max;
        List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors
        unsigned Num_Neighbors;                    // Number of neighbors found
        unsigned *Neighbor_IDs;                    // Array that holds final list of neighbors

        /* If we are near the edge of the cube then we need to adjust which
        particles we search through

        Note: Because unsigned integers rollover, we need to be careful to
        structure our tests such that they do not modify i j or k. For example,
        if k = 0 then check if k - Support_Radius < 0 will ALWAYS return
        false since 0 - Support_Radius = ~4 billion (rollover!). However,
        structuring the checks in this way makes them less readible, so I have
        included a logically equivalent (if rollover is ignored) if statement
        as a comment for each check */

        // i index (x coordinate) checks
        if(i < ((*this).Support_Radius)) { p_min = 0; }                                        // Same as if i - (*this).Support_Radius < 0
        else { p_min  = i - (*this).Support_Radius; }


        if(i > (X_SIDE_LENGTH - 1) - (*this).Support_Radius) { p_max = X_SIDE_LENGTH - 1; }  // Same as if(i + (*this).Support_Radius > X_SIDE_LENGTH -1)
        else { p_max = i + (*this).Support_Radius; }

        // j index (y coordinate) checks
        if(j < (*this).Support_Radius) { q_min = 0; }                                        // Same as if(j - (*this).Support_Radius < 0)
        else { q_min = j - (*this).Support_Radius; }

        if(j > (Y_SIDE_LENGTH - 1) - (*this).Support_Radius) { q_max = Y_SIDE_LENGTH - 1; }  // Same as if(j + (*this).Support_Radius > Y_SIDE_LENGTH - 1)
        else { q_max = j + (*this).Support_Radius; }

        // k index (z coordinate) checks
        if(k < (*this).Support_Radius) { r_min = 0; }                                        // Same as if(k - (*this).Support_Radius < 0)
        else { r_min = k - (*this).Support_Radius; }

        if(k > (Z_SIDE_LENGTH - 1) - (*this).Support_Radius) { r_max = Z_SIDE_LENGTH - 1; }  // Same as if(k + (*this).Support_Radius > Z_SIDE_LENGTH - 1)
        else { r_max = k + (*this).Support_Radius; }

        // Loop through potential neighbors, generate neighbor list
        for(p = p_min; p <= p_max; p++) {
          for(r = r_min; r <= r_max; r++) {
            for(q = q_min; q <= q_max; q++) {
              // a given particle is NOT its own neighbor
              if(i == p && j == q && k == r) { continue; }

              if(Are_Neighbors(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*(Y_SIDE_LENGTH) + j,
                               p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q))
                Particle_Neighbor_List.Add_Back(p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q);
            } // for(q = q_min; q <= q_max; q++) {
          } // for(r = r_min; r <= r_max; r++) {
        } // for(p = p_min; p <= p_max; p++) {

        /* Now that we have the neighbor list, we can make it into an array. To do
        this, we allocate an array whose length is equal to the length of the
        neighbor list. We then populate this array with the elements of the list
        and finally send this off to the particle (whose neighbors we found) */
        Num_Neighbors = Particle_Neighbor_List.Node_Count();
        Neighbor_IDs = new unsigned[Num_Neighbors];

        for(p = 0; p < Num_Neighbors; p++) { Neighbor_IDs[p] = Particle_Neighbor_List.Remove_Front(); }

        // Now sent the Neighbor list to the particle
        Set_Neighbors(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*(Y_SIDE_LENGTH) + j, Num_Neighbors, Neighbor_IDs);

        /* Now free Neighbor_IDs array for next particle! */
        delete [] Neighbor_IDs;
      } // for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {
} // void Body::Find_Neighbors_Cuboid(void) {



void Body::Remove_Neighbor(const unsigned i, const unsigned Remove_Neighbor_ID) {
  // This function is used to remove 1 neighbor from an existing particle.

  // To be able to remove a neighbor, the particle in question must have neighbors!
  if(Particles[i].Neighbors_Are_Set == false || Particles[i].Num_Neighbors == 0) {
    printf("Particle %d has no neighbors! We can't remove %d\n", i, Remove_Neighbor_ID);
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
  Neighbor arrays and point this particle's arrays to the new Neighbor
  arrays. */

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
    Vol_p = Particles[Particles[i].Neighbor_IDs[p]].Vol;                       //        : mm^3
    New_A += Dyadic_Product((Vol_p*New_Grad_W[p]), New_R[p]);                  // New shape tensor : unitless Tensor
  } // for(j = 0; j < Num_Neighbors; j++) {

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
