#if !defined(PARTICLE_DAMAGE)
#define PARTICLE_DAMAGE

#include "Particle_Helpers.h"

////////////////////////////////////////////////////////////////////////////////
// Damage methods

void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {
  printf("Particle %d is damaged.\n",P_In.ID);
  /* Find the set of particles that are in a sqrt(3) radius of the damaged
  particle, all of these should be removed. */
  unsigned int i,j,k;                                      // index variables
  const double ROOT_THREE = 1.73205081;                    // square root of 3           : inter-particle spacings
  const double Inter_Particle_Spacing = Particle::Inter_Particle_Spacing;      //        : mm

  unsigned int Num_Neighbors = P_In.Num_Neighbors;         // Neighbors of damaged particle
  unsigned int Neighbor_ID;                                // Current neighbor ID
  List Damaged_Neighbor_List;                              // A list of the ID's of the particles within a sqrt(3) radius of the damaged particle. These particles must go
  Damaged_Neighbor_List.Add_Front(P_In.ID);                     // Add Damaged particle to the lists
  unsigned int Num_Damaged_Neighbors = 1;                  // How many particles we have found within a sqrt(3) radius of the damaged particles


  Vector X = P_In.X;                                       // Reference position of damaged particle
  Vector X_Neighbor;

  /* Cycle through the damaged particle's neighbors, searching for ones
  inside of the sqrt(3) radius (that we intend to remove) */
  for(i = 0; i < Num_Neighbors; i++) {
    Neighbor_ID = P_In.Neighbor_IDs[i];
    X_Neighbor = Particles[i].X;

    /* if particle is within that sart(3) radius (in units of inter-particle
    of the damaged particle, then we need to remove it from the block. Add
    it to the Damaged particle list. Note that we add a little bit on to root
    three to account for roundoff errors/ensure that everything within a sqrt(3)
    radius is effected. */
    if(P_In.Mag_R[i] < ROOT_THREE*Inter_Particle_Spacing + .01) {
      Damaged_Neighbor_List.Add_Front((int)Neighbor_ID);
      Num_Damaged_Neighbors++;
    } // if(Magnitude(X_Neighbor - X) < ROOT_THREE) {
  } // for(i = 0; i < Num_Neighbors; i++) {

  // Copy damaged particle list to an array. Each damaged particle needs to be
  // designated as such (set each damaged particle's D to 1)
  unsigned int * Damaged_Particle_IDs = new unsigned int[Num_Damaged_Neighbors];
  for(i = 0; i < Num_Damaged_Neighbors; i++) {
    Damaged_Particle_IDs[i] = (unsigned int)Damaged_Neighbor_List.Remove_Front();
    Particles[Damaged_Particle_IDs[i]].D = 1;
  } // for(i = 0; i < Num_Damaged_Neighbors; i++) {

  /* Now that we know which particles need to be removed, we can causally remove
  them from the particles array. To do this, we need to find the set of all
  particles that have the damaged particle as a neighbor. Luckily, we know
  that for all particles A and B, if A is neighbor of B then B is a neighbor of
  A. We can begin with the damaged particle's Neighbor IDs. For each neighbor
  particle we can redo its neighbors list to exclude the damaged particle.
  Once we have don this, we can recalibrate the neighbor particle's members
  using the new reduced list. */
  unsigned int Damaged_Particle_ID;                        // ID of the particle we want to causally remove from the array
  unsigned int Damaged_Particle_Num_Neighbors;             // Number of neighbors of the particle that we're removing
  unsigned int * Damaged_Particle_Neighbors;               // Points to the neighbor list of the damaged particle

  // The Old and New variables corresond to the neighbor of a damaged particle
  unsigned int Old_Num_Neighbors;                          // for the neighbor of a damaged particle: Number of neighbors before removing damaged particle
  unsigned int * Old_Neighbors;                            // for the neighbor of a damaged particle: Old neighbor_IDs array

  unsigned int New_Num_Neighbors;                          // for the neighbor of a damaged particle: Number of neighbors now
  unsigned int * New_Neighbors;                            // for the neighbor of a damaged particle: New neighbor_IDs array
  unsigned int k_new;                                      // placement index for New_Neighbors (see explanation in for loop)

  // Cycle through the damaged particles
  for(i = 0; i < Num_Damaged_Neighbors; i++) {
    // For each damaged particle, get its neighbors
    Damaged_Particle_ID = Damaged_Particle_IDs[i];
    Damaged_Particle_Num_Neighbors = Particles[Damaged_Particle_ID].Num_Neighbors;
    Damaged_Particle_Neighbors = Particles[Damaged_Particle_ID].Neighbor_IDs;

    // Cycle through the neighbors of the damaged particle
    for(j = 0; j < Damaged_Particle_Num_Neighbors; j++) {
      /* For each neighbor of a damaged particle, remove the damaged particle
      from its Neighbor_IDs list then recalculate members like A_Inv, R, etc...
      (using set_neighbors). */

      Neighbor_ID = Damaged_Particle_Neighbors[j];
      Old_Num_Neighbors = Particles[Neighbor_ID].Num_Neighbors;
      Old_Neighbors = Particles[Neighbor_ID].Neighbor_IDs;

      New_Num_Neighbors = Old_Num_Neighbors-1;
      New_Neighbors = new unsigned int[New_Num_Neighbors];

      // Remove the damaged particle from neighboring particle's Neighbo_IDs array
      k_new = 0;
      for(k = 0; k < Old_Num_Neighbors; k++) {
        // Skip damaged particle
        if(Old_Neighbors[k] == Damaged_Particle_ID)
          continue;

        /* If a given neighbor is not damaged, we add it to the new neighbor
        list.

        Notice that we use the 'k_new' index in the new neighbor array and the
        'k' index of the old neighbor array. The reason is that we want the
        new array to skip the damaged element. before the damaged element, k
        and k_new are the same, after it k is one more than k_new. This way,
        by the end, k = Old_Num_Neighbors - 1 and k_new = New_Num_Neighbors - 1
        ( = Old_Num_Neighbors - 2)*/
        New_Neighbors[k_new] = Old_Neighbors[k];
        k_new++;
      } // for(k = 0; k < Num_Neighbors; k++) {

      //////////////////////////////////////////////////////////////////////////
      // Reset the neighbors for the 'neighboring' particle (reset its neighbor
      // members without the damaged particle)

      /* When we set new neighbors, the Set_Neighbors function will allocate
      new dynamic arrays for the Particle's members (for the W, Grad_W, etc..
      dynamic arrays). Thus, before we can set the new neighbors, we need to
      free the old dynamic arrays, thereby preventing a memory leak. */
      delete [] Particles[Neighbor_ID].R;                                      //        : mm
      delete [] Particles[Neighbor_ID].Mag_R;                                  //        : mm
      delete [] Particles[Neighbor_ID].W;                                      //        : unitless
      delete [] Particles[Neighbor_ID].Grad_W;                                 //        : mm^-1
      delete [] Particles[Neighbor_ID].Neighbor_IDs;

      // We need to set the 'Has_Neighbors' paramater to false. Otherwise, we
      // won't be able to set the neighbors.
      Particles[Neighbor_ID].Has_Neighbors = false;

      // Now we can reset the neighbors
      Particles[Neighbor_ID].Set_Neighbors(New_Num_Neighbors, New_Neighbors, Particles);

      // Now we can free the New_Neighbors array (it will be reallocated in the
      // next loop cycle, we free it to prevent a memory leak)
      delete [] New_Neighbors;
    } // for(j = 0; j < Damaged_Particle_Num_Neighbors; j++) {
  } // for(i = 0; i < Num_Damaged_Neighbors; i++) {

  /* Now that we've causally removed the damaged particles from the particles
  array, we need to make the damaged particle think it has no neighbors */
  for(i = 0; i < Num_Damaged_Neighbors; i++) {
    Particles[Damaged_Particle_IDs[i]].Num_Neighbors = 0;
  } // for(i = 0; i < Num_Damaged_Neighbors; i++) {

  // Free the 'Damaged_Particles_IDs' dynamic array to prevent memory leak
  delete [] Damaged_Particle_IDs;
} // void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {

#endif
