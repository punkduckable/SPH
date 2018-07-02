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
  /* Note: all tensor and vector members will be default initialized to zero
  (because that's how those class's constructors work) */
} // Particle::Particle(void) {

Particle::Particle(const Particle & P_In) {
  //printf("Particle copy constructor\n");

  /* Copy elements from P_In to local elements. Note that an object of one class
  is able to access private data of another object of that class.

  Note: because the Particle class contains pointers, we need to perform a deep
  copy of those pointers to ensure that the pointed too location isn't deleted
  when a temporary object is created and then destroyed by the copy constructor.
  */

  // Element wise copy of NON-POINTER members
  Vol = P_In.Vol;                                                              //        : mm^3
  Mass = P_In.Mass;                                                            //        : g

  X = P_In.X;                                                                  //        : mm
  x = P_In.x;                                                                  //        : mm
  vel = P_In.vel;                                                              //        : mm/s

  First_Iteration = P_In.First_Iteration;
  P = P_In.P;                                                                  //        : Mpa
  F = P_In.F;                                                                  //        : unitless

  /* Deep copy of pointer members.
  To do this, we need to give the new particle the same content as the origional
  neighbor list and Grad_W_Tilde arrays, but have these array's stored in a new
  memory location. This way, when the copy particle is deleted, it doesn't
  delete the origional particle's array.

  The new particle will have the same number of neighbors as the origional
  particle. The neighbor list and Grad_W_Tilde array's of the new particle
  should therefore be of length Num_Neighbors. Likewise, we need to copy the
  origional particle's array contents into the new particle's array conetents.
  We do this on an element by element basis. */
  Has_Neighbors = P_In.Has_Neighbors;
  Num_Neighbors = P_In.Num_Neighbors;

  Neighbor_IDs = new unsigned int[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  Mag_R = new double[Num_Neighbors];                                           //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1

  //Grad_W_Tilde = new Vector[Num_Neighbors];

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_IDs[j] = P_In.Neighbor_IDs[j];
    R[j] = P_In.R[j];                                                          //        : mm
    Mag_R[j] = P_In.Mag_R[j];                                                  //        : mm
    W[j] = P_In.W[j];                                                          //        : unitless
    Grad_W[j] = P_In.Grad_W[j];                                                //        : mm^-1
    //Grad_W_Tilde[j] = P_In.Grad_W_Tilde[j];
  } // for(unsinged int j = 0; j < Num_Neighbors; j++) {
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  //printf("Removing particle\n");

  // Note, we should only free the memory if it has been allocated.
  if(Has_Neighbors == true) {
    delete [] R;                                                               //        : mm
    delete [] Mag_R;                                                           //        : mm
    delete [] W;                                                               //        : unitless
    delete [] Grad_W;                                                          //        : mm^-1
    //delete [] Grad_W_Tilde;
    delete [] Neighbor_IDs;
  }
} // Particle::~Particle(void) {



////////////////////////////////////////////////////////////////////////////////
//  Particle equality

Particle & Particle::operator=(const Particle & P_In) {
  /* Copy elements from P_In to local elements. Note that an object of one class
  is able to access private data of another object of that class.

  Note: because the Particle class contains pointers, we need to perform a deep
  copy of those pointers to ensure that the pointed too location isn't deleted
  when a temporary object is created and then destroyed by the copy constructor.
  */

  // Element wise copy of NON-POINTER members
  Vol = P_In.Vol;                                                              //        : mm^3
  Mass = P_In.Mass;                                                            //        : g

  X = P_In.X;                                                                  //        : mm
  x = P_In.x;                                                                  //        : mm
  vel = P_In.vel;                                                              //        : mm/s

  First_Iteration = P_In.First_Iteration;
  P = P_In.P;                                                                  //        : Mpa
  F = P_In.F;                                                                  //        : unitless

  /* Deep copy of pointer members.
  To do this, we need to give the new particle the same content as the origional
  neighbor list and Grad_W_Tilde arrays, but have these array's stored in a new
  memory location. This way, when the copy particle is deleted, it doesn't
  delete the origional particle's array.

  The new particle will have the same number of neighbors as the origional
  particle. The neighbor list and Grad_W_Tilde array's of the new particle
  should therefore be of length Num_Neighbors. Likewise, we need to copy the
  origional particle's array contents into the new particle's array conetents.
  We do this on an element by element basis. */
  Has_Neighbors = P_In.Has_Neighbors;
  Num_Neighbors = P_In.Num_Neighbors;

  Neighbor_IDs = new unsigned int[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  Mag_R = new double[Num_Neighbors];                                           //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1
  //Grad_W_Tilde = new Vector[Num_Neighbors];

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_IDs[j] = P_In.Neighbor_IDs[j];
    R[j] = P_In.R[j];                                                          //        : mm
    Mag_R[j] = P_In.Mag_R[j];                                                  //        : mm
    W[j] = P_In.W[j];                                                          //        : unitless
    Grad_W[j] = P_In.Grad_W[j];                                                //        : mm^-1
    //Grad_W_Tilde[j] = P_In.Grad_W_Tilde[j];
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

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

  // Allocate memory for the Neighbor_IDs, and Grad_W_Tilde array
  Neighbor_IDs = new unsigned int[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  Mag_R = new double[Num_Neighbors];                                           //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1
  //Grad_W_Tilde = new Vector[Num_Neighbors];

  /* Now that we know our neighbors IDs, we can figure out everything that we
  want to know about them. We an set the Neighbor_IDs, r, R, W, and Grad_W
  members. These can be used to calculate the shape matrix (and its inverse)! */

  int Neighbor_ID;                               // Keep track of current particle
  double Vj;                                     // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,0,0,0,0,0,0};                   // Shape Tensor                         : unitless

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = Neighbor_ID_Array[j];          // Get Neighbor ID (index in Particles array)
    Neighbor_IDs[j] = Neighbor_ID;               // Set jth element of Neighbor_IDs member

    //r[j] = Particles[Neighbor_ID].x - x;       // Spacial displacement vector
    R[j] = Particles[Neighbor_ID].X - X;         // Reference displacement vector        : mm
    Vj = Particles[Neighbor_ID].Vol;             // Neighbor Volume                      : mm^3
    Mag_R[j] = R[j].Magnitude();                   // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amp*(h - Mag_R[j])*(h - Mag_R[j])*(h - Mag_R[j]);                             //        :unitless
    Grad_W[j] = -3*Shape_Function_Amp*((h - Mag_R[j])*(h - Mag_R[j]))*(R[j] / Mag_R[j]); //    : mm^-1

    // Add in the Current Neighbor's contribution to the Shape tensor
    A += Dyadic_Product((Vj*Grad_W[j]), R[j]);                                 //        : unitless
  } // for(unsigned int j = 0; j < N; j++) {

  // Now we can calculate A^(-1) from A.
  A_Inv = A^(-1);                                                              //        : unitless

  /* Now we can popuate the Grad_W_Tilde array
  for(unsigned int j = 0; j < N; j++) {
    Grad_W_Tilde[j] = A_Inv*Grad_W[j];
  } // for(unsigned int j = 0; j < N; j++) { */

  // Now that neighbors have been set, we set 'Has_Neighbors' to true
  Has_Neighbors = true;
} // void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_Id_List, const Particle *Particles) {

////////////////////////////////////////////////////////////////////////////////
// Friend functions (Update P, Update particle position)

void Update_P(Particle & P_In, const Particle * Particles, const double dt) {
  /* The purpose of this function is to calculate the First Piola-Kirchhoff
  stress tensor for the particle P_In.

  This function assumes that the position of each of P_In's neighbors has
  been updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its Grad_W_Tilde
  vector. Finally, it assumes that the static member variables k1, k2,
  and mu0 have been set.This function should not be called until these
  assumptions are valid.

  what are the arguments? This function accepts a Particle (P_In), a list of all
  particles in the current body, and the desired time step. This function uses
  these arguments to calculate P (the first Piola-Kirchhoff stress tensor) */

  /* First, let's set up the local variables that will be used to update the
  particle's position */
  double Vj;                                     // Volume of jth neighbor               : mm^3
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

  double Max_Principle_Stretch;
  const double Critical_Stretch = P_In.Critical_Stretch;

  Tensor A_Inv = P_In.A_Inv;                     // Inverse of shape tensor              : unitless
  const double Lame = P_In.Lame;                 // Lame paramater                       : Mpa
  const double mu0 = P_In.mu0;                   // Shear modulus                        : Mpa
  const double mu = P_In.mu;                     // Viscosity                            : Mpa*s
  const unsigned int Num_Neighbors = P_In.Num_Neighbors;

  Tensor F_Prime;                                // F time derivative                    : 1/s
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s
  Vector *Grad_W = P_In.Grad_W;                  // Pointer to P_In's Grad_W array.      : mm^-1
  Vector rj;                                     // Displacemtn vector of jth neighbor   : mm

  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) Vj Grad_W(Rj, h) */
  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = P_In.Neighbor_IDs[j];
    Vj = Particles[Neighbor_ID].Vol;                                           //        : mm^3
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm

    //F += Dyadic_Product(P_In.rj, Vj*P_In.Grad_W_Tilde[j]);
    F += Dyadic_Product(rj, Vj*Grad_W[j]);                                     //        : unitless
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

  // Deformation gradient with correction
  F *= A_Inv;                                                                  //        : unitless

  /* Calculate Damage. Here we calculate the 'damage'.
  to do this, we first need to find the principle stretch. To do this, we
  need to find the square root of the biggest eigenvalue of the (right)
  Cauchy Green strain tensor. Luckily, this tensor will be used for later
  calculations. */
  C = (F^(T))*F;                                 // Right Cauchy-Green strain tensor     : unitless
  J = Determinant(F);                            // J is det of F                        : unitless

  // Calculate current principle stretch.
  Max_Principle_Stretch = sqrt(Max_Eigenvalue(C,'F'));

  if(Max_Principle_Stretch > Critical_Stretch && Max_Principle_Stretch > P_In.Max_Stretch)
    P_In.Max_Stretch = Max_Principle_Stretch;


  /* Now that we have calculated the deformation gradient, we need to calculate
  the first Piola-Kirchhoff stess tensor. To do this, however, we need to
  find the Second Piola-Kirchhoff stress tensor and the Viscosity term. */


  S = mu0*I + (-mu0 + 2.*Lame*log(J))*(C^(-1));                                //        : Mpa

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
  P_In.P = (F*S + Visc)*A_Inv;                                                 //         : Mpa
  P_In.F = F;                                                                  //         : unitless

} // void Update_P(const Particle & P_In, const Particle * Particles, const double dt) {

void Update_Particle_Position(Particle & P_In, const Particle * Particles, const double dt) {
  /* This function assumes that every particle in the Particle's array has
  an updated P tensor. Likewise, it assumes that the E and alpha static
  member variables have been set. This function should not be run until
  these assumptions are valid. */

  Vector Force_Int{0,0,0};                       // Internal Force vector                : N
  Vector Force_Ext{0,0,0};                       // External/body force                  : N
  Vector Force_Hg{0,0,0};                        // Hour-glass force                     : N

  Vector acceleration;                           // acceleration vector                  : mm/s^2

  unsigned int Neighbor_ID;                      // ID of current neighbor particle (in paritlce's array)

  /* Jth particle variables */
  double Vj;                                     // Volume of jth particle               : mm^3
  Tensor P_j;                                    // First Piola-Kirchhoff stress tensor  : Mpa
  Tensor F_j;                                    // Deformation gradient                 : unitless
  Vector rj;                                     // Displacement vector                  : mm

  /* P_In aliases (ith particle variables).
  notice that P_In/P_i does not chnage throughout this function. Therefore,
  all P_In variables are declared as consts to avoid accidential modification */
  const double alpha = P_In.alpha;               // alpha static member                  : unitless
  const double E = P_In.E;                       // Hourglass stiffness                  : Mpa

  const double Vi = P_In.Vol;                    // Volume of P_In                       : mm^3
  const double Mass = P_In.Mass;                 // P_i's mass                           : g

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

    Vj = Particles[Neighbor_ID].Vol;                                           //        : mm^3
    P_j = Particles[Neighbor_ID].P;                                            //        : Mpa
    Force_Int += (Vi*Vj)*((P_i + P_j)*Grad_W[j]);                              //        : N

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

    F_j = Particles[Neighbor_ID].F;                                          //        : unitless
    delta_ji = Vector_Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj;//: mm

    /* Finally, we calculate the hour glass force. However, it should be
    noted that each term of Force_Hg is multiplied by -(1/2), E, alpha,
    and Vi. However, these four quantities are constants. We can therefore
    pull these multiplications out of the summations (thereby saving
    several thousand floating point operations per particle!)*/
    Force_Hg += (((Vj*W[j])/(Mag_R[j]*Mag_R[j]*Mag_rj))*
                (delta_ij + delta_ji))*(rj);
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {
  Force_Hg *= -.5*E*Vi*alpha;                                                  //        : N

  /* Compute acceleration of particle at new position a(t_i+1).
  Note that all the forces we have calculated have been in units of Newtons.
  Our mass is in units of grams and we want the acceleration in units of
  mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
  multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
  then we get acceleration in mm/s^2. */
  acceleration = ((1e+6)*(1./Mass))*(Force_Int + Force_Ext + Force_Hg);        //        : mm/s^2

  /* Now update the velocity, position vectors. This is done using the
  'leap-frog' integration scheme. However, during the first step of this
  scheme, we need to use forward euler to get the initial velocity.*/
  if(P_In.First_Iteration == true) {
    P_In.First_Iteration = false;
    P_In.vel += (dt/2.)*acceleration;            // velocity starts at t_i+1/2           : mm/s
  } //   if(P_In.First_Iteration == true) {

  P_In.x += dt*P_In.vel;                         // x_i+1 = x_i + dt*v_(i+1/2)           : mm
  P_In.vel += dt*acceleration;                   // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s

  P_In.Force_Int = Force_Int;
  P_In.Force_Ext = Force_Ext;
  P_In.Force_Hg = Force_Hg;
} // void Update_Particle_Position(Particle & P_In, const Particle * Particles, const double dt) {

bool Are_Neighbors(const Particle & P1, const Particle & P2) {
  /* This function checks if h > |Rj|. Here, Rj is simply the displacement of
  particle i relative to particle j: Rj = Xj - Xi. Xj = P1.X, Xi = P2.X. if
  h > |Rj| then P1 and P2 are in each other's support radius, so P1 is a
  neighbor of P2. */

  return ( P1.h > Magnitude(P1.X - P2.X));
}

////////////////////////////////////////////////////////////////////////////////
// Printing functions

void Particle::Print(void) const {
  // Print basic particle parameters.
  printf("Mass :   %f\n",Mass);
  printf("Volume: %f\n",Vol);
  printf("X:   ");
  X.Print();
  printf("x:   ");
  x.Print();
  printf("vel: ");
  vel.Print();
  printf("F:   \n");
  F.Print();
  printf("P:   \n");
  P.Print();
  printf("A^(-1)\n");
  A_Inv.Print();

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

    /* Print Grad_W_Tilde magnitudes
    printf("|Grad_W_Tilde| : {");
    for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
      printf("%5.3f, ",Magnitude(Grad_W_Tilde[i]));
    } // for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
    printf("%5.3f } \n", Magnitude(Grad_W_Tilde[Num_Neighbors-1])); // */

  } // if(Has_Neighbors == true) {
} // void Particle::Print(void) const {

void Print(const Particle & P_In) {
  P_In.Print();
} // void Print(const Particle & P_In) {

////////////////////////////////////////////////////////////////////////////////
// Generate particles

void Generate_Neighbor_Lists(const unsigned int Num_Particles, Particle * Particles) {
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
} // void Generate_Neighbor_Lists(const unsigned int Num_Particles, const Particle * Particles) {

void Generate_Neighbor_Lists_Box(const unsigned int Num_Particles, Particle * Particles,
                                 const unsigned int num_x, const unsigned int num_y, const unsigned int num_z,
                                 const unsigned int Support_Radius) {
  /* This function is a modified version of the Neighbor List generating
  function that is specialized for Box particle geometries. By box, I mean
  some kind of cuboid.

  Let us establish a few definitions:
  By a 'Layer' we mean a sheet of particles that all have the same x coordinate
  By a 'Vertical column' we mean a set of particles with the same x and y
  coordinates.
  By a 'Row' we mean a set of particles with the same x and z
  coordinates.

  This function assumes that the particles are stored in 'Vertical-Column'
  major, 'Layer' semi-major order. This menas that vertical columns of
  particles are stored in contiguous memory and that vertical columns in the
  same layer are stored in contiguous memory.

  Thus, if working with a cube with sidelength N, the (1,1,1) particle will be
  N*N particles away from the (2,1,1) partilce in the, N particles away from the
  (1,2,1) particle and 1 particle away from the (1,1,2) particle in the
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
  few extra arguments. the num_x, num_y, and num_z arguments specify the
  dimensions of the cuboid in the x, y, and z directions respectivly. Thus, if
  the cuboid has n layers, then num_x is n. If the cuboid has p particles in a
  vertical column then num_z is p. For a 100x50x200 cuboid of particles, num_x
  is 100, num_y is 50, and num_z is 200 */

  unsigned int i,j,k,p,q,r;            // Loop index variables
  unsigned int p_min, p_max, q_min, q_max, r_min, r_max;
  List Particle_Neighbor_List;         // Linked list to store known neighbors
  unsigned int Num_Neighbors;          // Number of neighbors found
  unsigned int *Neighbor_IDs;          // Array that holds final list of neighbors

  /* Cycle through the particles. For each particle, check if the particles near us in the cuboid grid
  are neighbors. */
  for(i = 0; i < num_x; i++) {
    for(j = 0; j < num_y; j++) {
      for(k = 0; k < num_z; k++) {

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
        if(i < Support_Radius)                   // Same as if(i - Support_Radius < 0).
          p_min = 0;
        else
          p_min  = i - Support_Radius;

        if(i > (num_x - 1) - Support_Radius)     // Same as if(i + Support_Radius > num_x -1)
          p_max = num_x - 1;
        else
          p_max = i + Support_Radius;

        // j index (y coordinate) checks
        if(j < Support_Radius)                   // Same as if(j - Support_Radius < 0)
          q_min = 0;
        else
          q_min = j - Support_Radius;

        if(j > (num_y - 1) - Support_Radius)     // Same as if(j + Support_Radius > num_y - 1)
          q_max = num_y - 1;
        else
          q_max = j + Support_Radius;

        // k index (z coordinate) checks
        if(k < Support_Radius)                   // Same as if(k - Support_Radius < 0)
          r_min = 0;
        else
          r_min = k - Support_Radius;

        if(k > (num_z - 1) - Support_Radius)     // Same as if(k + Support_Radius > num_z - 1)
          r_max = num_z - 1;
        else
          r_max = k + Support_Radius;

        // Loop through potential neighbors, generate neighbor list
        for(p = p_min; p <= p_max; p++) {
          for(q = q_min; q <= q_max; q++) {
            for(r = r_min; r <= r_max; r++) {
              // a given particle is NOT its own neighbor
              if(i == p && j == q && k == r)
                continue;

              if(Are_Neighbors(Particles[i*(num_z*num_y) + j*(num_z) + k] , Particles[p*(num_z*num_y) + q*(num_z) + r])) {
                Particle_Neighbor_List.Add_Back(p*(num_z*num_y) + q*(num_z) + r);
              }
            } // for(r = r_min; r <= r_max; r++) {
          } // for(q = q_min; q <= q_max; q++) {
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
        Particles[i*(num_z*num_y) + j*(num_z) + k].Set_Neighbors(Num_Neighbors, Neighbor_IDs, Particles);

        /* Now free Neighbor_IDs array for next particle! */
        delete [] Neighbor_IDs;
      } // for(k = 0; k < k_max; k++) {
    } // for(j = 0; j < j_max; j++) {
  } // for(i = 0; i < i_max; i++) {
} // void Generate_Neighbor_Lists_Box(const unsigned int Num_Particles, Particle * Particles,
  //                                  const unsigned int num_x, const unsigned int num_y, const unsigned int num_z,
  //                                  const unsigned int Support_Radius) {

#endif
