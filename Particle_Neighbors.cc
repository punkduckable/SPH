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

  return ( P1.h > Magnitude(P1.X - P2.X));
} // bool Particle_Helpers::Are_Neighbors(const Particle & P1, const Particle & P2) {

void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles, Particle * Particles) {
  unsigned int i,j;                    // Loop index variables
  List Particle_Neighbor_List;         // Linked list to store known neighbors
  unsigned int Num_Neighbors;          // Number of neighbors found
  unsigned int *Neighbor_IDs;          // Array that holds final list of neighbors

  // Cycle through the particles
  for(i = 0; i < Num_Particles; i++) {

    /* For each particle, cycle through the potential neighbors (every particle) */
    for(j = 0; j < Num_Particles; j++) {
      // ith particle is not its own neighbor.
      if(j == i)
        continue;

      // Test if jth particle is inside support radius of ith particle
      if(Are_Neighbors(Particles[i], Particles[j])) {
        Particle_Neighbor_List.Add_Back(j);
      } // if(Are_Neighbors(Particles[i], Particles[j])) {
    } // for(unsigned int j = 0; j < Num_Particles; j++) {

    /* Now that we have the neighbor list, we can make it into an array. To do
    this, we allocate an array whose length is equal to the length of the
    neighbor list. We then populate this array with the elements of the list
    and finally send this off to the particle (whose neighbors we found) */
    Num_Neighbors = Particle_Neighbor_List.Node_Count();
    Neighbor_IDs = new unsigned int[Num_Neighbors];

    for(j = 0; j < Num_Neighbors; j++) {
      Neighbor_IDs[j] = Particle_Neighbor_List.Remove_Front();
    } // for(j = 0; j < Num_Neighbors; j++) {

    // Now sent the Neighbor list to the particle
    Particles[i].Set_Neighbors(Num_Neighbors, Neighbor_IDs, Particles);

    /* Now free Neighbor_IDs array for next particle! */
    delete [] Neighbor_IDs;
  } // for(unsigned int i = 0; i < Num_Particles; i++) {
} // void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles, Particle * Particles) {

void Particle_Helpers::Find_Neighbors_Box(Particle & P_In, Particle * Particles) {
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

  unsigned int i = P_In.ijk[0], j = P_In.ijk[1], k = P_In.ijk[2];
  unsigned int p,q,r;                  // Loop index variables
  unsigned int p_min, p_max, q_min, q_max, r_min, r_max;
  List Particle_Neighbor_List;         // Linked list to store known neighbors
  unsigned int Num_Neighbors;          // Number of neighbors found
  unsigned int *Neighbor_IDs;          // Array that holds final list of neighbors

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
          Particle_Neighbor_List.Add_Back(p*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + r*(Y_SIDE_LENGTH) + q);
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

#endif
