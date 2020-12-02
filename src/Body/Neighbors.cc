#include "Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Array.h"
#include "Errors.h"
#include <assert.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////////////
// Neighbor methods!

void Body::Set_Neighbors(const unsigned i, const Array<unsigned> & Neighbor_IDs_In) {
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



void Body::Find_Neighbors(void) {
  unsigned i,j;                              // Loop index variables
  List<unsigned> Particle_Neighbor_List;     // Linked list to store known neighbors

  // Cycle through the particles
  for(i = 0; i < Num_Particles; i++) {

    /* For each particle, cycle through the potential neighbors (every particle) */
    for(j = 0; j < Num_Particles; j++) {
      // ith particle is not its own neighbor.
      if(j == i) { continue; }

      // Test if jth particle is inside support radius of ith particle. If so,
      // add P_j to P_i's neighbor list.
      if(Are_Neighbors(i, j)) { Particle_Neighbor_List.Push_Back(j); }
    } // for(unsigned j = 0; j < Num_Particles; j++) {

    /* Now that we have the neighbor ID list, we can make it into an array.
    This is done using the Array class' list constructor. See Array.h */
    Array<unsigned> Neighbor_IDs(Particle_Neighbor_List);

    // Now sent the Neighbor list to the particle
    Set_Neighbors(i, Neighbor_IDs);
  } // for(unsigned i = 0; i < Num_Particles; i++) {
} // void Body::Find_Neighbors(void) {



void Body::Find_Neighbors_Box(void) {
  if((*this).Is_Box == false) {
    char Buf[500];
    sprintf(Buf,
            "Not A Box Exception: thrown by Body::Find_Neighbors_Box\n"
            "Body %s tried to use this function, but %s is not a box! This function\n"
            "can only be called by boxes!\n",
            (*this).Name.c_str(), (*this).Name.c_str());
    throw Not_A_Box(Buf);
  } // if((*this).Is_Box == false) {

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

  const unsigned d_max = floor((*this).Support_Radius / (*this).Inter_Particle_Spacing);

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
  // This function is used to remove 1 neighbor from an existing particle.

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
