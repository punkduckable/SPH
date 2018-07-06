#if !defined(PARTICLE_SOURCE)
#define PARTICLE_SOURCE

#include "Particle.h"
#include "Tensor.h"
#include "Vector.h"
#include "List.h"

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  //printf("Particle default constructor \n");
  Has_Neighbors = false;
  Num_Neighbors = 0;
  Vol = 0;                                                                     //        : mm^3
  Mass = 0;                                                                    //        : g

  // Now randomly set critical stress
  unsigned seed = std::rand();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(1.3,.05);
  Stretch_Critical = distribution(generator);
} // Particle::Particle(void) {

Particle::Particle(const Particle & P_In) {
  // Do nothing copy constructor
  printf("Bad! You tried using the particle copy constructor. Don't do that!\n");
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  //printf("Removing particle\n");

  // Note, we should only free the memory if it has been allocated.
  if(Has_Neighbors == true) {
    delete [] R;                                                               //        : mm
    delete [] Mag_R;                                                           //        : mm
    delete [] W;                                                               //        : unitless
    delete [] Grad_W;                                                          //        : mm^-1
    delete [] Neighbor_IDs;
  }
} // Particle::~Particle(void) {



////////////////////////////////////////////////////////////////////////////////
//  Particle equality

Particle & Particle::operator=(const Particle & P_In) {
  // Do nothing = operator overload.
  printf("You can't use equality with particles!\n");
  return *this;
}



////////////////////////////////////////////////////////////////////////////////
// Particle setup methods:
// Set Neighbors

void Particle::Set_Neighbors(const unsigned int N, const unsigned int * Neighbor_ID_Array, const Particle * Particles) {
  /* First check if this particle already has neighbors. This function should
  only be called if the neighbors have not been set. The reason for this is
  that this method allocates pointers. If the pointers have already been set,
  then allocating them again will cause a memory leak.  */
  if(Has_Neighbors == true) {
    printf("Neigbors already set, returning.\n");
    return;
  }

  /* Next, check that N > 0. if N = 0, then there are no neighbors to add. */
  if(N == 0) {
    printf("You didn't supply any neighbors!\n");
    return;
  }

  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Dynamic arrays
  Neighbor_IDs = new unsigned int[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  Mag_R = new double[Num_Neighbors];                                           //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1

  /* Now that we know our neighbors IDs, we can figure out everything that we
  want to know about them. We an set the Neighbor_IDs, r, R, W, and Grad_W
  members. These can be used to calculate the shape matrix (and its inverse)! */

  int Neighbor_ID;                               // Keep track of current particle
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,
           0,0,0,
           0,0,0};                               // Shape Tensor (zero initialized)      : unitless

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = Neighbor_ID_Array[j];          // Get Neighbor ID (index in Particles array)
    Neighbor_IDs[j] = Neighbor_ID;               // Set jth element of Neighbor_IDs member

    // Calculate displacement vectors
    R[j] = Particles[Neighbor_ID].X - X;         // Reference displacement vector        : mm
    Mag_R[j] = R[j].Magnitude();                 // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amp*(h - Mag_R[j])*(h - Mag_R[j])*(h - Mag_R[j]);                             //        :unitless
    Grad_W[j] = -3*Shape_Function_Amp*((h - Mag_R[j])*(h - Mag_R[j]))*(R[j] / Mag_R[j]); //    : mm^-1

    // Add in the Current Neighbor's contribution to the Shape tensor
    V_j = Particles[Neighbor_ID].Vol;            // Neighbor Volume                      : mm^3
    A += Dyadic_Product((V_j*Grad_W[j]), R[j]);                                //        : unitless
  } // for(unsigned int j = 0; j < N; j++) {

  // Now we can calculate A^(-1) from A.
  A_Inv = A^(-1);                                                              //        : unitless

  // Now that neighbors have been set, we set 'Has_Neighbors' to true
  Has_Neighbors = true;
} // void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_Id_List, const Particle *Particles) {

////////////////////////////////////////////////////////////////////////////////
// Friend functions (Update P, Update particle position)

void Update_P(Particle & P_In, Particle * Particles, const double dt) {
  // Check if particle is damaged (if so, we skip this particle)
  if(P_In.D >= 1)
    return;

  /* The purpose of this function is to calculate the First Piola-Kirchhoff
  stress tensor for the particle P_In.

  This function assumes that the position of each of P_In's neighbors has
  been updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its dynamic array
  members. Finally, it assumes that the static member variables k1, k2,
  and mu0 have been set.This function should not be called until these
  assumptions are valid.

  what are the arguments? This function accepts a Particle (P_In), a list of all
  particles in the current body, and the desired time step. This function uses
  these arguments to calculate P (the first Piola-Kirchhoff stress tensor) */

  /* First, let's set up the local variables that will be used to update the
  particle's position */
  double V_j;                                     // Volume of jth neighbor               : mm^3
  unsigned int Neighbor_ID;                      // Index of jth neighbor.

  Tensor F = {0,0,0,
              0,0,0,
              0,0,0};                            // Deformation gradient                 : unitless

  Tensor C;                                      // Richt-Cauchy stress tensor           : unitless
  double J;                                      // Deformation gradient determinant     : unitless
  Tensor S;                                      // Second Poila-Kirchhoff stress tensor : Mpa
  Tensor I = {1,0,0,
              0,1,0,
              0,0,1};                            // Identity tensor

  double Stretch_Max_Principle;
  const double Tau = Particle::Tau;

  const double Lame = Particle::Lame;            // Lame paramater                       : Mpa
  const double mu0 = Particle::mu0;              // Shear modulus                        : Mpa
  const double mu = Particle::mu;                // Viscosity                            : Mpa*s
  const unsigned int Num_Neighbors = P_In.Num_Neighbors;

  Tensor F_Prime;                                // F time derivative                    : 1/s
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s
  Vector *Grad_W = P_In.Grad_W;                  // Pointer to P_In's Grad_W array.      : mm^-1
  Vector rj;                                     // Displacemtn vector of jth neighbor   : mm

  //////////////////////////////////////////////////////////////////////////////
  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) V_j Grad_W(Rj, h) */
  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = P_In.Neighbor_IDs[j];
    V_j = Particles[Neighbor_ID].Vol;                                          //        : mm^3
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm

    //F += Dyadic_Product(P_In.rj, V_j*P_In.Grad_W[j])*A^(-1);
    F += Dyadic_Product(rj, V_j*Grad_W[j]);                                    //        : unitless
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

  // Deformation gradient with correction
  F *= P_In.A_Inv;                                                             //        : unitless

  //////////////////////////////////////////////////////////////////////////////
  /* Calculate Damage:

  To do this, we first need to find the principle stretch. To do this, we
  need to find the square root of the biggest eigenvalue of the (right)
  Cauchy Green strain tensor. Luckily, this tensor will be used for later
  calculations. */
  C = (F^(T))*F;                                 // Right Cauchy-Green strain tensor     : unitless
  J = Determinant(F);                            // J is det of F                        : unitless

  // Calculate current principle stretch
  Stretch_Max_Principle = sqrt(Max_Eigenvalue(C,'F'));

  // If this stretch is greater than max stretch, update particle's Max stretch.
  P_In.Stretch_M = Stretch_Max_Principle;
  if(Stretch_Max_Principle > P_In.Stretch_H)
    P_In.Stretch_H = Stretch_Max_Principle;

  // if Max is greater than crticial and the particle is in the rip zone then
  // start adding damage
  if(P_In.ijk[1] == Y_SIDE_LENGTH/2 || P_In.ijk[1] == Y_SIDE_LENGTH/2-1 || P_In.ijk[1] == Y_SIDE_LENGTH/2+1) {
    if(P_In.Stretch_H > P_In.Stretch_Critical)
      P_In.D = exp(((P_In.Stretch_H - P_In.Stretch_Critical)*(P_In.Stretch_H - P_In.Stretch_Critical))/(Tau*Tau)) - 1;

    // If particle is fully damaged, remove it from array.
    if(P_In.D >= 1) {
      Remove_Damaged_Particle(P_In, Particles);
      return;
    } // if(P_In.D >= 1) {
  } // if(P_In.ijk[1] == Y_SIDE_LENGTH/2 || P_In.ijk[1] == Y_SIDE_LENGTH/2-1) {

  //////////////////////////////////////////////////////////////////////////////
  /* Now that we have calculated the deformation gradient, we need to calculate
  the first Piola-Kirchhoff stess tensor. To do this, however, we need to
  find the Second Piola-Kirchhoff stress tensor and the Viscosity term. */


  S = (1-P_In.D)*(mu0*I + (-mu0 + 2.*Lame*log(J))*(C^(-1)));                        //        : Mpa

  /* Calculate viscosity tensor:
  To do this, we need to calculate the deformation gradient. Luckily, at this
  point, the F member of P_In (the deformation gradient that P_In contains) has
  not been updated. Thus, P_In's F member is the old strain tensor. Therefore,
  we can use the F that we calculated above as the 'new' deformation tensor,
  F(t), and P_In.F as the 'old' deformation tensor, F(t-dt). We can then use
  the forward difference approximation of the derivative to get an approximation
  for F_Prine. */
  F_Prime = (1./dt)*(F - P_In.F);                                              //        : s^-1
  L = F_Prime*(F^(-1));                                                        //        : s^-1
  Visc = (J*mu)*(L + (L^(T))*(F^(-T)));                                        //        : Mpa

  /* Calculate P (First Piola-Kirchhoff stress tensor), send it and F to P_In */
  P_In.P = (F*S + Visc)*P_In.A_Inv;                                            //         : Mpa
  P_In.F = F;                                                                  //         : unitless
  //P_In.Visc = Visc*P_In.A_Inv;                                                 // For debugging

} // void Update_P(const Particle & P_In, Particle * Particles, const double dt) {

void Update_x(Particle & P_In, const Particle * Particles, const double dt) {
  // Check if particle is damaged (if so, we skip this particle)
  if( P_In.D >= 1)
    return;

  /* This function assumes that every particle in the Particle's array has
  an updated P tensor. Likewise, it assumes that the E and alpha static
  member variables have been set. This function should not be run until
  these assumptions are valid. */

  P_In.Force_Int = {0,0,0};                      // Internal Force vector                : N
  P_In.Force_Ext = {0,0,0};                      // External/body force                  : N
  P_In.Force_Hg = {0,0,0};                       // Hour-glass force                     : N
  //P_In.Force_Visc = {0,0,0};                     // For debugging

  Vector acceleration;                           // acceleration vector                  : mm/s^2

  unsigned int Neighbor_ID;                      // ID of current neighbor particle (in paritlce's array)

  /* Jth particle variables */
  double V_j;                                    // Volume of jth particle               : mm^3
  Tensor P_j;                                    // First Piola-Kirchhoff stress tensor  : Mpa
  Tensor F_j;                                    // Deformation gradient                 : unitless
  Vector rj;                                     // Displacement vector                  : mm

  /* P_In aliases (ith particle variables).
  notice that P_In/P_i does not chnage throughout this function. Therefore,
  all P_In variables are declared as consts to avoid accidential modification */
  const double alpha = Particle::alpha;          // alpha static member                  : unitless
  const double E = Particle::E;                  // Hourglass stiffness                  : Mpa

  const double V_i = P_In.Vol;                   // Volume of P_In                       : mm^3

  const Tensor P_i = P_In.P;                     // First Piola-Kirchhoff stress tensor  : Mpa
  const Tensor F_i = P_In.F;                     // Deformation gradient                 : unitless

  const unsigned int Num_Neighbors = P_In.Num_Neighbors; // Number of neighbors of P_In
  const Vector * R = P_In.R;                     // Reference displacement array         : mm
  const double * Mag_R = P_In.Mag_R;             // Mag of reference displacement array  : mm
  const double * W = P_In.W;                     // Shape function array                 : unitless
  const Vector * Grad_W = P_In.Grad_W;           // Grad_W array                         : 1/mm

  /* Hour glass variables */
  double Mag_rj;                                                               //        : mm
  double delta_ij;                                                             //        : mm
  double delta_ji;                                                             //        : mm

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    // Update Neighbor
    Neighbor_ID = P_In.Neighbor_IDs[j];

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate Internal force */

    /* Note, each term in the internal force sum is multiplied by Vi. If  we
    we were to multiply through by Vi in this loop, we'd peerform this operation
    Num_Neighbors times. By moving it out of the summation (mutiplying Force_Int
    by Vi after the loop) we reduce the number of multiplications to 1, thereby
    reducing the number of FLOPs required to calculate the internal force and
    speeding up the program. */

    V_j = Particles[Neighbor_ID].Vol;                                          //        : mm^3
    P_j = Particles[Neighbor_ID].P;                                            //        : Mpa
    P_In.Force_Int += (V_j)*((P_i + P_j)*Grad_W[j]);                           //        : N
    //P_In.Force_Visc += (V_j)*((P_In.Visc + Particles[Neighbor_ID].Visc)*Grad_W[j]); // For debugging

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate external Force */

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate Hour Glass force */

    /* Here we calculate delta_ij.
    Before discussing this, let us establish the following definitions
          r_ij = rj - ri
          R_ij = Rj - Ri
          rk = spacial position of kth particle
          Rk = ref position of kth particle
          Fk = deformation gradient of kth particle)
    delta_ij is given by,
          delta_ij = (Error_ij dot r_ij)/|r_ij|
    where
          Error_ij = Fi*R_ij - r_ij.
    since the dot product is distributive, notice that
          delta_ij = ((Fi*R_ij - r_ij) dot (r_ij))/|r_ij|
                  = (Fi*R_ij dot r_ij)/|r_ij| - (r_ij dot r_ij)/|r_ij|
                  = (Fi*R_ij dot r_ij)/|r_ij| - |r_ij|^2/|r_ij|
                  = (Fi*R_ij dot r_ij)/|r_ij| - |r_ij|
    Calcualating delta this way actually uses fewer floating point operations
    and should therefore perform better. */
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm
    Mag_rj = Magnitude(rj);                                                    //        : mm
    delta_ij = Vector_Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj;             //        : mm

    /* Here we calculate delta_ji.
          delta_ji = ( Error_ji dot r_ji )/|r_ji|
    With
          Error_ji = F_j*R_ji - r_ji
    Notice that we need the jth particles deformation gradient to calculate
    Error_ji. We also need the jth particles R_ji and r_ij. We could calculate all
    of these, but we can save some time by making a few clever observations.
    From the definion of R,
          R_ji = X_i - X_j = -(X_j - X_i) = -R_ij
    Likewise, from the defintion of r, we can decduce
          r_ij = -r_ji
    Thus,
        Error_ji = F_j*R_ji - r_ji
                 = F_j*(-R_ij) + r_ij
                 = -F_j*(R_ij) + r_ij
    Then, using the fact that the dot product is distributive, we get
          delta_ji = ( Error_ji dot r_ji ) / |r_ji|
                   = ( Error_ji dot -r_ij ) / |r_ij|
                   = -( Error_ji dot r_ij ) / |r_ij|
                   = -( (-F_j*R_ij + r_ij ) dot r_ij )/|r_ij|
                   = -( -F_j*R_ij dot r_ij) / |r_ij| - (r_ij dot r_ij)/|r_ij|
                   = (F_j*R_ij dot r_ij) / |r_ij| - |r_ij|^2/|r_ij|
                   = (F_j*R_ij dot r_ij) / |r_ij| - |r_ij|
    Computing delta_ji this way uses fewer arithmetic operations and should
    therefore improve performnace.
    */

    F_j = Particles[Neighbor_ID].F;                                            //        : unitless
    delta_ji = Vector_Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj;//: mm

    /* Finally, we calculate the hour glass force. However, it should be
    noted that each term of Force_Hg is multiplied by -(1/2), E, alpha,
    and Vi. However, these four quantities are constants. We can therefore
    pull these multiplications out of the summations (thereby saving
    several thousand floating point operations per particle!)*/
    P_In.Force_Hg += (((V_j*W[j])/(Mag_R[j]*Mag_R[j]*Mag_rj))*
                (delta_ij + delta_ji))*(rj);
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {
  P_In.Force_Hg *= -.5*E*V_i*alpha;    // Each term in F_Hg is multiplied by this. Pulling it out of sum improved runtime
  P_In.Force_Int *= V_i;               // Each term in F_Int sum is multiplied by Vi, pulling it out of sum improved runtime
  //P_In.Force_Visc *= V_i;              // for debugging

  /* Compute acceleration of particle at new position a(t_i+1).
  Note that all the forces we have calculated have been in units of Newtons.
  Our mass is in units of grams and we want the acceleration in units of
  mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
  multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
  then we get acceleration in mm/s^2. */
  acceleration = ((1e+6)*(1./P_In.Mass))*(P_In.Force_Int +
                                          P_In.Force_Ext +
                                          P_In.Force_Hg);                      //        : mm/s^2

  /* Now update the velocity, position vectors. This is done using the
  'leap-frog' integration scheme. However, during the first step of this
  scheme, we need to use forward euler to get the initial velocity.*/
  if(P_In.First_Time_Step == true) {
    P_In.First_Time_Step = false;
    P_In.vel += (dt/2.)*acceleration;            // velocity starts at t_i+1/2           : mm/s
  } //   if(P_In.First_Time_Step == true) {

  P_In.x += dt*P_In.vel;                         // x_i+1 = x_i + dt*v_(i+1/2)           : mm
  P_In.vel += dt*acceleration;                   // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s
} // void Update_x(Particle & P_In, const Particle * Particles, const double dt) {

////////////////////////////////////////////////////////////////////////////////
// Neighbor methods!

bool Are_Neighbors(const Particle & P1, const Particle & P2) {
  /* This function checks if h > |Rj|. Here, Rj is simply the displacement of
  particle i relative to particle j: Rj = Xj - Xi. Xj = P1.X, Xi = P2.X. if
  h > |Rj| then P1 and P2 are in each other's support radius, so P1 is a
  neighbor of P2. */

  return ( P1.h > Magnitude(P1.X - P2.X));
} // bool Are_Neighbors(const Particle & P1, const Particle & P2) {

void Find_Neighbors(const unsigned int Num_Particles, Particle * Particles) {
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
} // void Find_Neighbors(const unsigned int Num_Particles, const Particle * Particles) {

void Find_Neighbors_Box(Particle & P_In, Particle * Particles) {
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
} // void Find_Neighbors(Particle & P_In, Particle * Particles) {



////////////////////////////////////////////////////////////////////////////////
// Damage methods

void Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {
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
} // void Remove_Damaged_Particle(Particle & P_In, Particle * Particles) {

////////////////////////////////////////////////////////////////////////////////
// Printing functions

void Particle::Print(void) const {
  // Print basic particle parameters.
  printf("X:   ");
  X.Print();
  //printf("x:   ");
  //x.Print();
  //printf("vel: ");
  //vel.Print();
  //printf("F:   \n");
  //F.Print();
  //printf("P:   \n");
  //P.Print();
  //printf("A^(-1)\n");
  //A_Inv.Print();
  printf("Num Neighbors: %d\n",Num_Neighbors);
  printf("F_Int = ");
  Force_Int.Print();
  printf("F_Visc = ");
  //Force_Visc.Print();                          // For debugging
  //printf("F_Hg = ");
  Force_Hg.Print();
  printf("\n");

  // If we have neighbors, print neighbor information
  if(Has_Neighbors == true) {
    //unsigned int i;                    // Loop index variable

    /* Print neighbor ID's
    printf("Neighbor ID's  : {");
    for(i = 0; i < Num_Neighbors-1; i++) {
      printf("%5d, ",Neighbor_IDs[i]);
    } // for(i = 0; i < Num_Neighbors-1; i++) {
    printf("%5d } \n", Neighbor_IDs[Num_Neighbors-1]);  // */

    /* Print Grad_W magnitudes
    printf("%p\n",Grad_W);
    printf("|Grad_W|       : {");
    for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
      printf("%5.3f, ",Magnitude(Grad_W[i]));
    } // for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
    printf("%5.3f } \n", Magnitude(Grad_W[Num_Neighbors-1])); // */

  } // if(Has_Neighbors == true) {
} // void Particle::Print(void) const {

void Print(const Particle & P_In) {
  P_In.Print();
} // void Print(const Particle & P_In) {

#endif
