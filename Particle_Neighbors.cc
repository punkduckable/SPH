#if !defined(PARTICLE_NEIGHBORS)
#define PARTICLE_NEIGHBORS

#include "Particle_Helpers.h"

////////////////////////////////////////////////////////////////////////////////
// Neighbor methods!

bool Particle_Helpers::Are_Neighbors(const Particle & P1, const Particle & P2) {
  /* This function checks if h > |Rj|. Here, Rj is simply the displacement of
  particle i relative to particle j: Rj = Xj - Xi. Xj = P1.X, Xi = P2.X. if
  h > |Rj| then P1 and P2 are in each other's support radius, so P1 is a
  neighbor of P2. */

  return ( Particle::h > Magnitude(P1.Get_X() - P2.Get_X()));
} // bool Particle_Helpers::Are_Neighbors(const Particle & P1, const Particle & P2) {

void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles, Particle * Particles) {
  unsigned int i,j;                              // Loop index variables
  List<unsigned int> Particle_Neighbor_List;     // Linked list to store known neighbors
  unsigned int Num_Neighbors;                    // Number of neighbors found
  unsigned int *Neighbor_IDs;                    // Array that holds final list of neighbors

  // Cycle through the particles
  for(i = 0; i < Num_Particles; i++) {

    /* For each particle, cycle through the potential neighbors (every particle) */
    for(j = 0; j < Num_Particles; j++) {
      // ith particle is not its own neighbor.
      if(j == i)
        continue;

      // Test if jth particle is inside support radius of ith particle. If so,
      // add P_j to P_i's neighbor list.
      if(Are_Neighbors(Particles[i], Particles[j]))
        Particle_Neighbor_List.Add_Back(Particles[j].Get_ID());
    } // for(unsigned int j = 0; j < Num_Particles; j++) {

    /* Now that we have the neighbor list, we can make it into an array. To do
    this, we allocate an array whose length is equal to the length of the
    neighbor list. We then populate this array with the elements of the list
    and finally send this off to the particle (whose neighbors we found) */
    Num_Neighbors = Particle_Neighbor_List.Node_Count();
    Neighbor_IDs = new unsigned int[Num_Neighbors];

    for(j = 0; j < Num_Neighbors; j++)
      Neighbor_IDs[j] = Particle_Neighbor_List.Remove_Front();

    // Now sent the Neighbor list to the particle
    Particles[i].Set_Neighbors(Num_Neighbors, Neighbor_IDs, Particles);

    /* Now free Neighbor_IDs array for next particle! */
    delete [] Neighbor_IDs;
  } // for(unsigned int i = 0; i < Num_Particles; i++) {
} // void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles, Particle * Particles) {

void Particle_Helpers::Find_Neighbors_Box(Particle & P_In, Particle * Particles, const unsigned X_SIDE_LENGTH, const unsigned Y_SIDE_LENGTH, const unsigned Z_SIDE_LENGTH) {
  /* This function is a modified version of the Neighbor List generating
  function that is specialized for Box particle geometries. By box, I mean
  some kind of cuboid.

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
  This function is used just like the Generage_Neighbor_List function but with a
  few extra arguments. the X_SIDE_LENGTH, Y_SIDE_LENGTH, and Z_SIDE_LENGTH arguments specify the
  dimensions of the cuboid in the x, y, and z directions respectivly. Thus, if
  the cuboid has n layers, then X_SIDE_LENGTH is n. If the cuboid has p particles in a
  vertical column then Z_SIDE_LENGTH is p. For a 100x50x200 cuboid of particles, X_SIDE_LENGTH
  is 100, Y_SIDE_LENGTH is 50, and Z_SIDE_LENGTH is 200 */

  unsigned int i = P_In.Get_i(), j = P_In.Get_j(), k = P_In.Get_k();
  unsigned int p,q,r;                             // Loop index variables
  unsigned int p_min, p_max, q_min, q_max, r_min, r_max;
  List<unsigned int> Particle_Neighbor_List;     // Linked list to store known neighbors
  unsigned int Num_Neighbors;                    // Number of neighbors found
  unsigned int *Neighbor_IDs;                    // Array that holds final list of neighbors
  const unsigned int SUPPORT_RADIUS = Particle::Support_Radius;

  /* If we are near the edge of the cube then we need to adjust which
  particles we search through

  Note: Because unsigned integers rollover, we need to be careful to
  structure our tests such that they do not modify i j or k. For example,
  if k = 0 then check if k - SUPPORT_RADIUS < 0 will ALWAYS return
  false since 0 - SUPPORT_RADIUS = ~4 billion (rollover!). However,
  structuring the checks in this way makes them less readible, so I have
  included a logically equivalent (if rollover is ignored) if statement
  as a comment for each check */

  // i index (x coordinate) checks
  if(i < SUPPORT_RADIUS)                         // Same as if(i - SUPPORT_RADIUS < 0).
    p_min = 0;
  else
    p_min  = i - SUPPORT_RADIUS;

  if(i > (X_SIDE_LENGTH - 1) - SUPPORT_RADIUS)   // Same as if(i + SUPPORT_RADIUS > X_SIDE_LENGTH -1)
    p_max = X_SIDE_LENGTH - 1;
  else
    p_max = i + SUPPORT_RADIUS;

  // j index (y coordinate) checks
  if(j < SUPPORT_RADIUS)                         // Same as if(j - SUPPORT_RADIUS < 0)
    q_min = 0;
  else
    q_min = j - SUPPORT_RADIUS;

  if(j > (Y_SIDE_LENGTH - 1) - SUPPORT_RADIUS)   // Same as if(j + SUPPORT_RADIUS > Y_SIDE_LENGTH - 1)
    q_max = Y_SIDE_LENGTH - 1;
  else
    q_max = j + SUPPORT_RADIUS;

  // k index (z coordinate) checks
  if(k < SUPPORT_RADIUS)                         // Same as if(k - SUPPORT_RADIUS < 0)
    r_min = 0;
  else
    r_min = k - SUPPORT_RADIUS;

  if(k > (Z_SIDE_LENGTH - 1) - SUPPORT_RADIUS)   // Same as if(k + SUPPORT_RADIUS > Z_SIDE_LENGTH - 1)
    r_max = Z_SIDE_LENGTH - 1;
  else
    r_max = k + SUPPORT_RADIUS;

  // Loop through potential neighbors, generate neighbor list
  for(p = p_min; p <= p_max; p++) {
    for(r = r_min; r <= r_max; r++) {
      for(q = q_min; q <= q_max; q++) {
        // a given particle is NOT its own neighbor
        if(i == p && j == q && k == r)
          continue;

        if(Are_Neighbors(Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*(Y_SIDE_LENGTH) + j], Particles[p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q]))
          Particle_Neighbor_List.Add_Back(Particles[p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q].Get_ID());
      } // for(q = q_min; q <= q_max; q++) {
    } // for(r = r_min; r <= r_max; r++) {
  } // for(p = p_min; p <= p_max; p++) {

  /* Now that we have the neighbor list, we can make it into an array. To do
  this, we allocate an array whose length is equal to the length of the
  neighbor list. We then populate this array with the elements of the list
  and finally send this off to the particle (whose neighbors we found) */
  Num_Neighbors = Particle_Neighbor_List.Node_Count();
  Neighbor_IDs = new unsigned int[Num_Neighbors];

  for(p = 0; p < Num_Neighbors; p++) {
    Neighbor_IDs[p] = Particle_Neighbor_List.Remove_Front();
  } // for(j = 0; j < Num_Neighbors; j++) {

  // Now sent the Neighbor list to the particle
  P_In.Set_Neighbors(Num_Neighbors, Neighbor_IDs, Particles);

  /* Now free Neighbor_IDs array for next particle! */
  delete [] Neighbor_IDs;
} // void Particle_Helpers::Find_Neighbors_Box(Particle & P_In, Particle * Particles) {


void Particle_Helpers::Remove_Neighbor(Particle & P_In, const unsigned int Remove_Neighbor_ID, const Particle * Particles) {
  // This function is used to remove 1 neighbor from an existing particle.

  // To be able to remove a neighbor, we need to have neighbors!
  if(P_In.Neighbors_Are_Set == false || P_In.Num_Neighbors == 0)
    printf("Particle %d has no neighbors! We can't remove %d\n", P_In.ID, Remove_Neighbor_ID);

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

  unsigned int i, j = -1;                              // index variables
  double V_j;

  // New neighbor arrays
  unsigned int *New_Neighbor_IDs = new unsigned int[P_In.Num_Neighbors - 1];   //        : unitless
  Vector *New_R = new Vector[P_In.Num_Neighbors - 1];                          //        : mm Vector
  double *New_Mag_R = new double[P_In.Num_Neighbors - 1];                      //        : mm
  double *New_W = new double[P_In.Num_Neighbors - 1];                          //        : 1/(mm^3)
  Vector *New_Grad_W = new Vector[P_In.Num_Neighbors - 1];                     //        : 1/(mm^4) Vector
  Tensor New_A{0,0,0,
               0,0,0,
               0,0,0};

  for(i = 0; i < P_In.Num_Neighbors; i++) {
    // Check if ith neighbor ID matches Remove_Neighbor_ID
    if(P_In.Neighbor_IDs[i] == Remove_Neighbor_ID)
      continue;

    // If not, then this is a new neighbor, increment j.
    j++;

    // Check if j == Num_Neighbors - 1. If this is the case, then Remove_Particle_ID
    // was NOT one of this particle's neighbors!
    if(j == P_In.Num_Neighbors - 1) {
      printf("%d was not a neighbor of %d\n",Remove_Neighbor_ID, P_In.ID);
      return;
    }

    // Copy old Neighbor data to New Neighbor arrays.
    New_Neighbor_IDs[j] = P_In.Neighbor_IDs[i];                                //        : unitless
    New_R[j] = P_In.R[i];                                                      //        : mm Vector
    New_Mag_R[j] = P_In.Mag_R[i];                                              //        : mm
    New_W[j] = P_In.W[i];                                                      //        : 1/(mm^3)
    New_Grad_W[j] = P_In.Grad_W[i];                                            //        : 1/(mm^4) Vector

    // Calculate New shape tensor.
    V_j = Particles[P_In.Neighbor_IDs[j]].Vol;                                 //        : mm^3
    New_A += Dyadic_Product((V_j*New_Grad_W[j]), New_R[j]);          // New shape tensor : unitless Tensor
  } // for(i = 0; i < Num_Neighbors; i++) {

  // Now that we have our new neighbor arrays, we can replace/delete the old
  // neighbor arrays
  delete [] P_In.Neighbor_IDs;
  delete [] P_In.R;
  delete [] P_In.Mag_R;
  delete [] P_In.W;
  delete [] P_In.Grad_W;

  P_In.Neighbor_IDs = New_Neighbor_IDs;                                        //        : unitless
  P_In.R = New_R;                                                              //        : mm Vector
  P_In.Mag_R = New_Mag_R;                                                      //        : mm
  P_In.W = New_W;                                                              //        : 1/(mm^3)
  P_In.Grad_W = New_Grad_W;                                                    //        : 1/(mm^4) Vector

  // Decrement number of neighbors by 1.
  P_In.Num_Neighbors--;

  // Now we can calculate the new A^(-1) from New_A.
  P_In.A_Inv = New_A^(-1);                                                     //        : unitless Tensor

} // void Particle_Helpers::Remove_Neighbor(Particle & P_In, const unsigned int Remove_Neighbor_ID, const Particle * Particles, const unsigned X_SIDE_LENGTH, const unsigned Y_SIDE_LENGTH, const unsigned Z_SIDE_LENGTH) {

#endif
