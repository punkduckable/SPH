#include "Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Array.h"

////////////////////////////////////////////////////////////////////////////////
// Damage methods

void Body::Remove_Damaged_Particles(List<unsigned> & Damaged_Particle_List) {
  while(Damaged_Particle_List.Get_Num_Nodes() != 0) {

  unsigned p = Damaged_Particle_List.Pop_Back();

  #if defined(DAMAGE_MONITOR)
    printf("Particle %d is damaged. ", p);
    Particles[p].X.Print();
  #endif

  /* Particle p (which we will refer to as P) is damaged. We want to causally
  remove P from the Particles array. This means that the position of P can no
  longer impact the stress or motion of any particles in the array. How do we do
  this? By removing P from every other particle's neighbor lists! But we must do
  more. Let P_i, P_j and P be particles with the following reference
  configuration

         P_i              P                P_j
         ()               ()               ()

  Suppose that P is damaged. Further suppose that P_i and P_j are neighbors.
  Since P is damaged, if the body breaks at P, then there can't be any
  connection from P_i to P_j. Therefore, if P is infront of P_j, then P_i
  shouldn't be neighbors with P_j. We say that P_j is in P's 'shadow region'
  Basically, if P_i were a light and P were some obstical, then any particle
  P_j in the 'shadow region' would be covered by the shadow of P. More
  More formally, for any two particles, P_i and P_j, if the ray between P_i and
  P_j passes through the damaged particle, then P_j is in the 'shadow region' of
  P. For any particle P_i, we remove all neighbors that are damaged or in the
  shadow region of a damaged particle.

  Throughout this discission, all rays (position vectors) are with respect to
  REFERENCE positions. We do not care about spatial positions.

  We can actually greatly speed up this process by making a fundamental
  observation: If the ray between P_i and P_j passes through P, then P
  must be closer to P_i than P_j. This means that P must ALSO be a neighbor
  of P_i! Since neighborship is symmetric, this means that P_i must be one of
  P's neighbors. Therefore, for any damaged particle P, the only particles
  that can be in P's shadow region are neighbors of P. Thus, we don't
  really need to search through the entire Particle's array; rather, we only
  need to search through the damaged particle's neighbor list!

  How do we actually detect if the ray from P_i to P_j intersects the damaged
  particle? Each particle has a 'radius' r. We treat the particle as a sphere
  centered at the Particle's reference position with radius r. This radius is
  typically half the inter-particle spacing. If the Ray from P_i to P_j
  interescts the sphere, then P_j is in P's shadow region and is therefore
  no longer one of P_i's neighbors.

  So how do we detect if the ray intersects the sphere? There are three
  necessairry conditions for ray-sphere intersection:

  For the fist condition, Let R_i R_j and R_In denote the reference position
  vectors of P_i P_j and P respectivly. If the ray from R_i to R_j intersects
  P's sphere, then the scalar projection of (R_In-R_i) onto (R_j - R_i) must
  be positive. Why is this? Well, if this projection is negaitve then
  (R_In - R_i) must go in away from (R_j - R_i) and therefore can not intersect
  the sphere. Therefore, our first condition for intersection is:

         (R_In - R_i) dot (R_j - R_i) > 0


  For the second condition, we need to check the lengths of (R_In - R_i) and
  (R_j - R_i). If the ray (R_j - R_i) is to pass through P's sphere, then
  (R_j - R_i) must be longer than (R_In - R_i)? To understand why, consider the
  following diagrams:


                                 Intersection:
             R_i                       ____                      R_j
             ()----------------------/-----\---------------------()
                                    |  ()  |
                                    \_____/
                                      R_In


                               No Intersection:

             R_i            R_j        ____
             ()-------------()       /     \
                                    |  ()  |
                                    \_____/
                                      R_In

  Even if (R_i - R_j) is pointed in the same direction as (R_In - R_i) (the
  scalar projection of (R_In - R_i) onto (R_i - R_j) is positive), intersection
  occurs only if Notice that if (R_i - R_j) is LONGER than (R_i - R_In).
  Therefore, our second necessairry condition for intersction is:

         |R_In - R_i| < |R_j - R_i|

  For the third conidition: Assuming that the two conditions avbove check out
  then we need to find the shortest distance between a point on (R-i - R_j) and
  R_In. Consider the following diagram:

  R_i                                                                  R_j
  ()___________________________________________________________________()
    -----_____                         L|
              -----_____                | d
                        -----_____      |
                                  -----()
                                       R_In

  Notice that the rays (R_i - R_j), d, and (R_In - R_i) form a right triangle.
  If the ray (R_j - R_i) intersects P's sphere then d^2 < r^2 where r is the
  radius of P's sphere. However, notice that

         d^2 = ((R_In - R_i) dot (R_In - R_i)) - ((R_In - R_i) dot (R_i - R_j)/|R_i - R_j|))^2

  With this, we can check our third and final condition,

         d^2 < r^2.

  If all three conditions are passed, then P_j is in the shadow region of P;
  thus, we remove P_i and P_j's neighbor status. */

  const double r_Squared = Particles[p].Radius*Particles[p].Radius;            //        : mm^2

  // Particle i (P_j) paramaters
  unsigned Pi_ID;                            // ID of P_i
  List<unsigned> Pi_New_Neighbor_List;       // List of all of P_i's neighbors that are not damaged or in Particles[p]'s shadow region

  // Particle j (P_j) paramaters
  unsigned Pj_ID;                            // ID of P_j

  // Rays (Vectors between particle's Reference positions)
  Vector Rj_Ri;                                  // R_j - R_i                            : mm Vector
  Vector RIn_Ri;                                 // R_In - R_i                           : mm Vector
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
  passes through Particles[i]'s' shadow region. If so, we remove its neighbor
  status. This is done by redoing P_i's neighbors list to exclude the damaged
  particle. Once we have done this, we can recalibrate the neighbor particle's
  members using the new reduced list. */
  unsigned Num_Neighbors = Particles[p].Num_Neighbors;
  for(unsigned i = 0; i < Num_Neighbors; i++) {
    /* For each neighbor of the damaged particle, P_i, we need to remove the
    neirhbors of Particles[i] as well as any of P_i's neighbors that are in the
    shadow region from P_i's neighbor list. Thus, we need to cycle through P_i's
    neighbors. */
    Pi_ID = Particles[p].Neighbor_IDs[i];

    RIn_Ri = Particles[p].X - Particles[Pi_ID].X;                                      //        : mm Vector

    for(unsigned j = 0; j < Particles[Pi_ID].Num_Neighbors; j++) {
      Pj_ID = Particles[Pi_ID].Neighbor_IDs[j];

      // First, check if P_j is the damaged particle
      if(Pj_ID == Particles[p].ID) { continue; }


      //////////////////////////////////////////////////////////////////////////
      // Checks: Now we check if P_j is in Particles[i]'s shadow region
      Rj_Ri = Particles[Pj_ID].X - Particles[Pi_ID].X;                         //        : mm Vecor


      //////////////////////////////////////////////////////////////////////////
      // First check:
      // Check if the scalar projection of R_In - R_i onto Rj - Rj is positive
      // or negative.

      if( Dot_Product(RIn_Ri, Rj_Ri) <= 0 ) {
        // If not positive, then Ri and Rj must be neighbors, add P_j to P_i's
        // new neighbor list.
        Pi_New_Neighbor_List.Push_Front(Pj_ID);
        continue;
      } // if( Dot_Product(RIn_Ri, Rj_Ri) <= 0 ) {


      //////////////////////////////////////////////////////////////////////////
      // Second Check:
      // If RIn_Ri dot Rj_Ri is positive, then we need to perform the second
      // check. To do this, we make sure that the magnitude of RIn_Rj is less
      // than the magnitude of Rj_Ri.

      RIn_Ri_Dot_RIn_Ri = Dot_Product(RIn_Ri, RIn_Ri);                         //        : mm
      Rj_Ri_Dot_Rj_Ri = Dot_Product(Rj_Ri, Rj_Ri);                             //        : mm

      if(RIn_Ri_Dot_RIn_Ri >= Rj_Ri_Dot_Rj_Ri) {
        // if |RIn_Ri| > |Rj_Ri| then P_j is NOT in the shadow region of
        // Particles[i]. thus, P_j is still a neighbor of P_i.
        Pi_New_Neighbor_List.Push_Front(Pj_ID);
        continue;
      } // if(RIn_Ri_Dot_RIn_Ri >= Rj_Ri_Dot_Rj_Ri) {


      //////////////////////////////////////////////////////////////////////////
      // Third check:
      // If the first two checks are passed, then we can run the third check
      // by calculating d^2 and comparing it to r^2.

      // Calculate d^2
      RIn_Ri_Dot_Rj_Ri = Dot_Product(RIn_Ri, Rj_Ri);                           //        : mm

      d_Squared = Dot_Product(RIn_Ri, RIn_Ri) - (RIn_Ri_Dot_Rj_Ri * RIn_Ri_Dot_Rj_Ri) / Rj_Ri_Dot_Rj_Ri;        // mm

      // Check if d^2 < r^2
      if(d_Squared <= r_Squared) {
        /* if so, then P_j is in P_i's shadow region! This means that P_i is
        no longer a neighbor of P_j, and that P_j is no longer a neighbor of
        P_i. To things need to happen for this to work, P_i needs to stop
        being a neighbor with P_j and P_j needs to stop being a neighbor
        with P_i.

        To ensure that P_j is not one of P_i's neighbors, we simply do NOT
        include P_j on P_i's new neighbor's list.

        To ensure that P_i is not one of P_j's neighbors, we need to remove
        P_i from P_j's neighbor list. This is done by the 'Remove Neighbor'
        member function.

        Removing this neighbor may give Pj a non-invertible A matrix. If this
        is the case, then calling Remove_Neighbor will throw a Singular_Matrix
        exception. If this happens, then we need to damage Pj and add it to
        the damaged particle list. */
        try { (*this).Remove_Neighbor(Pj_ID, Pi_ID); }
        catch(Singular_Matrix & Er_In) {
          Particles[Pj_ID].Set_D(1);
          Damaged_Particle_List.Push_Front(Pj_ID);
        } // catch(Singular_Matrix & Er_In) {

        continue;
      } // if(d_squared < r_squared) {

      else {
        // otherwise, P_i and P_j are neighbors.
        Pi_New_Neighbor_List.Push_Front(Pj_ID);
        continue;
      } // else {
    } // for(j = 0; j < Particles[Pi_ID].Num_Neighbors; j++) {


    ////////////////////////////////////////////////////////////////////////////
    /* We now have a complete new neighbor list for P_i (with Particles[i] and the
    shadow region particles removed). We can now begin the process of redoing
    P_i's member variables. */

    // To begin, convert the New Neirghbor List to an Array
    Array<unsigned> Pi_New_Neighbors(Pi_New_Neighbor_List);


    /* When we set new neighbors, the Set_Neighbors function will allocate
    new dynamic arrays for the Particle's members (for the W, Grad_W, etc..
    dynamic arrays). Thus, before we can set the new neighbors, we need to
    free the old dynamic arrays, thereby preventing a memory leak. */
    delete [] Particles[Pi_ID].R;                                              //        : mm Vector
    delete [] Particles[Pi_ID].Mag_R;                                          //        : mm
    delete [] Particles[Pi_ID].W;                                              //        : 1/(mm^3)
    delete [] Particles[Pi_ID].Grad_W;                                         //        : 1/(mm^4) Vector
    delete [] Particles[Pi_ID].Neighbor_IDs;

    // We need to set the 'Neighbors_Are_Set' paramater to false. Otherwise, we
    // won't be able to reset P_i's neighbors
    Particles[Pi_ID].Neighbors_Are_Set = false;

    /* Now we can reset the neighbors. This may throw a Singular_Matrix
    exception. If this happens, then we need to damage Pi and it to the Damaged
    particle list */
    try { Set_Neighbors(Pi_ID, Pi_New_Neighbors); }
    catch(Singular_Matrix & Er_In) {
      Particles[Pi_ID].Set_D(1);
      Damaged_Particle_List.Push_Front(Pi_ID);
    } // catch(Singular_Matrix & Er_In) {
  } // for(i = 0; i < Num_Neighbors; i++) {

  /* Now that we've causally removed the damaged particle from the particles
  array, we need to make that it think that it has no neighbors */
  Particles[p].Num_Neighbors = 0;
  } // while(Damaged_Particle_List.Get_Num_Nodes() != 0) {
} // void Body::Remove_Damaged_Particles(List<unsigned> & Damaged_Particle_List) {
