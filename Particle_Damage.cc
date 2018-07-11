#if !defined(PARTICLE_DAMAGE)
#define PARTICLE_DAMAGE

#include "Particle_Helpers.h"

////////////////////////////////////////////////////////////////////////////////
// Damage methods

void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {
  printf("Particle %d is damaged.\n",P_In.ID);
  /* Particle P_In is damaged. This basically means that we must causally remove
  it from the Particles array. How do we do this? By modifying neighbor lists!
  A damaged particle creates a 'shadow region'. For any two particles, P_i and
  P_j, if the ray between P_i and P_j passes through the damaged particle, then
  P_j is considered to be in the 'shawdow region' of P_i. For any particle P_i,
  if any of its neighbors are damaged or are in the shadow region of a damaged
  particle then we remove those neighbors.

  Throughout this discission, all rays (position vectors) are with respect to
  REFERENCE positions. We do not care about the current positions.

  We can actually greatly speed up this process by making a fundamental
  observation: If the ray between P_i and P_j passes through P_In, then P_In
  must be closer to P_i than P_j. This means that P_In must ALSO be an neighbor
  of P_i! Since neighborship is symmetric, this means that P_i must be one of
  P_In's neighbors. Thus, we don't really need to search through the entire
  Particle's array; rather, we only need to search through the damaged particle's
  neighbor list!

  How do we actually detect if the ray from P_i to P_j intersects the damaged
  particle? Each particle has a 'radius'. We treat the particle as the sphere
  enclosed by that radius (that is centered at the particle's position). The
  radius is usually the inter-particle-spacing. if the Ray from P_i to P_j
  interescts the sphere, then we remove P_i and P_j's neighbor status.

  So how do we detect if the ray intersects the sphere? Let R_i R_j and R_In
  denote the position vectors of P_i P_j and P_In respectivly. if the ray between
  R_i and R_j intersects the sphere around R_In, then the scalar projection
  of (R_In-R_i) onto (R_j - R_i) must be positive. If its negative then the ray
  between P_i and P_j goes away from P_In, meaning that the ray can not
  intersect P_In's sphere. Therefore, our first condition for intersection is:

        (R_In - R_i) dot (R_j - R_i) > 0


  Assuming that the above condition checks out, then we at least know that the
  ray is pointed in the correct direction. Further, consider the following
  diagram:

  R_i                                                                  R_j
  ()___________________________________________________________________()
    -----_____                         L|
              -----_____                | d
                        -----_____      |
                                  -----()
                                       R_In

  Notice that the rays (R_i - R_j), d, and (R_In - R_i) form a right triangle.
  If the ray (R_j - R_i) intersects P_In's sphere then d^2 < r^2 where r is the
  radius of P_In's sphere. However, notice that
        d^2 = ((R_In - R_i) dot (R_In - R_i)) - ((R_In - R_i) dot (R_i - R_j)/|R_i - R_j|))^2
  With this, we can check the second condition, d^2 < r^2. If true, then P_j is
  in the shadow region of P_In; thus, we remove P_i and P_j's neighbor status.
  */

  unsigned int i,j,k;                           // index variables
  const double r_Squared = P_In.Inter_Particle_Spacing*P_In.Inter_Particle_Spacing/4.;      //        : mm^2

  // Damaged particle (P_In) paramaters
  unsigned int Damaged_Particle_ID = P_In.ID;    // ID of the damaged particle

  // Particle i (P_j) paramaters
  unsigned int Pi_ID;                            // ID of P_i
  unsigned int Pi_New_Num_Neighbors;             // Number of neighbors of P_i
  List Pi_New_Neighbor_List;                     // List of all of P_i's neighbors that are not damaged or in P_In's shadow region
  unsigned int * Pi_New_Neighbors;               // Stores P_i's new neighbors.

  // Particle j (P_j) paramaters
  unsigned int Pj_ID;                            // ID of P_j

  // Rays (Vectors between particle's Reference positions)
  Vector Rj_Ri;                                  // R_j - R_i                  : mm
  Vector RIn_Ri;                                 // R_In - R_i                 : mm
  double RIn_Ri_Dot_RIn_Ri;                      // Dot product of RIn_Ri and RIn_Ri (= |RIn_Ri|^2). Used to determine intersection
  double RIn_Ri_Dot_Rj_Ri;                       // Dot product of RIn_Ri and Rj_Ri. Used to calculate d^2.            : mm^2
  double Rj_Ri_Dot_Rj_Ri;                        // Dot product of Rj_Ri and Rj_Ri (= |Rj-Ri|^2). Used to calculate d^2: mm^2

  // d (orthogonal distance from ray to sphere)
  double d_Squared;                                                            //        : mm^2

  /* Here we causally remove the damaged particle from the Particles array.
  To do this, we need to find the set of all particles that are neighbors with
  the damaged particle. Luckily, since Neighborship is symmetric, for all
  particles A and B, if A is neighbor of B then B is a neighbor of A. Thus, we
  simply cycle through each if the damaged particle's neighbors. For each
  neighbor, P_i we check if the ray between P_i and each of its neighbors, P_j
  passes through P_In shadow region. If so, we remove its neighbor status. This is
  done by redoing P_i's neighbors list to exclude the damaged particle.
  Once we have done this, we can recalibrate the neighbor particle's members
  using the new reduced list. */
  for(i = 0; i < P_In.Num_Neighbors; i++) {
    /* For each neighbor of the damaged particle, P_i, we need to remove both P_In
    and any of P_i's neighbors that are in the shadow region from P_i's neighbor
    list. Thus, we need to cycle through P_i's neighbors. */
    Pi_ID = P_In.Neighbor_IDs[i];

    RIn_Ri = P_In.X - Particles[Pi_ID].X;

    for(j = 0; j < Particles[Pi_ID].Num_Neighbors; j++) {
      Pj_ID = Particles[Pi_ID].Neighbor_IDs[j];

      // First, check if P_j is the damaged particle
      if(Pj_ID == Damaged_Particle_ID)
        continue;

      //////////////////////////////////////////////////////////////////////////
      // Checks: Now we check if P_j is in P_In's shadow region
      Rj_Ri = Particles[Pj_ID].X - Particles[Pi_ID].X;

      //////////////////////////////////////////////////////////////////////////
      // First check:
      // Check if the scalar projection of R_In - R_i onto Rj - Rj is positive
      // or negative.

      if( Vector_Dot_Product(RIn_Ri, Rj_Ri) <= 0) {
        // If not positive, then Ri and Rj must be neighbors, add P_j to P_i's
        // new neighbor list.
        Pi_New_Neighbor_List.Add_Front(Pj_ID);
        continue;
      } // if( Vector_Dot_Product(RIn_Ri, Rj_Ri) <= 0) {

      //////////////////////////////////////////////////////////////////////////
      // Second Check:
      // If RIn_Ri dot Rj_Ri is positive, then we need to perform the second
      // check. To do this, we make sure that the magnitude of RIn_Rj is less
      // than the magnitude of Rj_Ri.

      RIn_Ri_Dot_RIn_Ri = Vector_Dot_Product(RIn_Ri, RIn_Ri);
      Rj_Ri_Dot_Rj_Ri = Vector_Dot_Product(Rj_Ri, Rj_Ri);

      if(RIn_Ri_Dot_RIn_Ri >= Rj_Ri_Dot_Rj_Ri) {
        // if |RIn_Ri| > |Rj_Ri| then P_j is NOT in the shadow region of P_In.
        // thus, P_j is still a neighbor of P_i.
        Pi_New_Neighbor_List.Add_Front(Pj_ID);
        continue;
      } // if(RIn_Ri_Dot_RIn_Ri >= Rj_Ri_Dot_Rj_Ri) {

      //////////////////////////////////////////////////////////////////////////
      // Third check:
      // If the first two checks are passed, then we can run the third check
      // by calculating d^2 and comparing it to r^2.

      // Calculate d^2
      RIn_Ri_Dot_Rj_Ri = Vector_Dot_Product(RIn_Ri, Rj_Ri);

      d_Squared = Vector_Dot_Product(RIn_Ri, RIn_Ri) - (RIn_Ri_Dot_Rj_Ri * RIn_Ri_Dot_Rj_Ri) / Rj_Ri_Dot_Rj_Ri;

      // Check if d^2 < r^2
      if(d_Squared <= r_Squared) {
        // if so, then P_j is in P_i's shadow region!
        continue;
      } // if(d_squared < r_squared) {

      else {
        // otherwise, P_i and P_j are neighbors.
        Pi_New_Neighbor_List.Add_Front(Pj_ID);
        continue;
      } // else {
    } // for(j = 0; j < Particles[Pi_ID].Num_Neighbors; j++) {

    ////////////////////////////////////////////////////////////////////////////
    /* We now have a complete new neighbor list for P_i (with P_In and the
    shadow region particles removed). We can now begin the process of redoing
    P_i's member variables. */
    Pi_New_Num_Neighbors = Pi_New_Neighbor_List.Node_Count();
    Pi_New_Neighbors = new unsigned int[Pi_New_Num_Neighbors];

    for(k = 0; k < Pi_New_Num_Neighbors; k++) {
        Pi_New_Neighbors[k] = Pi_New_Neighbor_List.Remove_Front();
    } // for(k = 0; k < Pi_New_Num_Neighbors; k++) {

    /* When we set new neighbors, the Set_Neighbors function will allocate
    new dynamic arrays for the Particle's members (for the W, Grad_W, etc..
    dynamic arrays). Thus, before we can set the new neighbors, we need to
    free the old dynamic arrays, thereby preventing a memory leak. */
    delete [] Particles[Pi_ID].R;                                              //        : mm
    delete [] Particles[Pi_ID].Mag_R;                                          //        : mm
    delete [] Particles[Pi_ID].W;                                              //        : unitless
    delete [] Particles[Pi_ID].Grad_W;                                         //        : mm^-1
    delete [] Particles[Pi_ID].Neighbor_IDs;

    // We need to set the 'Has_Neighbors' paramater to false. Otherwise, we
    // won't be able to set the neighbors.
    Particles[Pi_ID].Has_Neighbors = false;

    // Now we can reset the neighbors
    Particles[Pi_ID].Set_Neighbors(Pi_New_Num_Neighbors, Pi_New_Neighbors, Particles);

    // Now we can free the New_Neighbors array (it will be reallocated in the
    // next loop cycle, we free it to prevent a memory leak)
    delete [] Pi_New_Neighbors;

  } // for(i = 0; i < P_In.Num_Neighbors; i++) {

  /* Now that we've causally removed the damaged particle from the particles
  array, we need to make that it think that it has no neighbors */
  P_In.Num_Neighbors = 0;
} // void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {

#endif
