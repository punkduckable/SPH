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
  //r = new Vector[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1
  //Grad_W_Tilde = new Vector[Num_Neighbors];

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_IDs[j] = P_In.Neighbor_IDs[j];
    //r[j] = P_In.r[j]
    R[j] = P_In.R[j];                                                          //        : mm
    W[j] = P_In.W[j];                                                          //        : unitless
    Grad_W[j] = P_In.Grad_W[j];                                                //        : mm^-1
    //Grad_W_Tilde[j] = P_In.Grad_W_Tilde[j];
  } // for(unsinged int j = 0; j < Num_Neighbors; j++) {
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  //printf("Removing particle\n");

  // Note, we should only free the memory if it has been allocated.
  if(Has_Neighbors == true) {
    //delete [] r;
    delete [] R;                                                               //        : mm
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
  //r = new Vector[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1
  //Grad_W_Tilde = new Vector[Num_Neighbors];

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_IDs[j] = P_In.Neighbor_IDs[j];
    //r[j] = P_In.r[j]
    R[j] = P_In.R[j];                                                          //        : mm
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
  //r = new Vector[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1
  //Grad_W_Tilde = new Vector[Num_Neighbors];

  /* Now that we know our neighbors IDs, we can figure out everything that we
  want to know about them. We an set the Neighbor_IDs, r, R, W, and Grad_W
  members. These can be used to calculate the shape matrix (and its inverse)! */

  int Neighbor_ID;                               // Keep track of current particle
  double Mag_Rj;                                 // magnitude of Rj vector               : mm
  double Vj;                                     // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,0,0,0,0,0,0};                   // Shape Tensor                         : unitless

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = Neighbor_ID_Array[j];          // Get Neighbor ID (index in Particles array)
    Neighbor_IDs[j] = Neighbor_ID;               // Set jth element of Neighbor_IDs member

    //r[j] = Particles[Neighbor_ID].x - x;       // Spacial displacement vector
    R[j] = Particles[Neighbor_ID].X - X;         // Reference displacement vector        : mm
    Vj = Particles[Neighbor_ID].Vol;             // Neighbor Volume                      : mm^3
    Mag_Rj = R[j].Magnitude();                   // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amp*(h - Mag_Rj)*(h - Mag_Rj)*(h - Mag_Rj);                             //        :unitless
    Grad_W[j] = -3*Shape_Function_Amp*((h - Mag_Rj)*(h - Mag_Rj))*(R[j] / Mag_Rj); //    : mm^-1

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
  int Neighbor_ID;                               // Index of jth neighbor.

  Tensor F = {0,0,0,0,0,0,0,0,0};                // Deformation gradient                 : unitless
  Tensor A_Inv = P_In.A_Inv;                     // Inverse of shape tensor              : unitless

  Tensor C;                                      // Richt-Cauchy stress tensor           : unitless
  Tensor S;                                      // Second Poila-Kirchhoff stress tensor : Mpa
  Tensor I = {1,0,0,
              0,1,0,
              0,0,1};                            // Identity tensor

  double Lame = P_In.Lame;                       // Lame paramater                       : Mpa
  double mu0 = P_In.mu0;                         // Shear modulus                        : Mpa
  double mu = P_In.mu;                           // Viscosity                            : Mpa*s
  double  J;                                     // Deformation gradient determinant     : unitless

  Tensor F_Prime;                                // F time derivative                    : 1/s
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s
  Vector *Grad_W = P_In.Grad_W;                  // Pointer to P_In's Grad_W array.      : mm^-1
  Vector rj;                                     // Displacemtn vector of jth neighbor   : mm

  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) Vj Grad_W(Rj, h) */
  for(unsigned int j = 0; j < P_In.Num_Neighbors; j++) {
    Neighbor_ID = P_In.Neighbor_IDs[j];
    Vj = Particles[Neighbor_ID].Vol;                                           //        : mm^3
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm

    //F += Dyadic_Product(P_In.rj, Vj*P_In.Grad_W_Tilde[j]);
    F += Dyadic_Product(rj, Vj*Grad_W[j]);                                     //        : unitless
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

  // Deformation gradient with correction
  F = F*A_Inv;                                                                 //        : unitless

  /* Now that we have calculated the deformation gradient, we need to calculate
  the first Piola-Kirchhoff stess tensor. To do this, however, we need to
  find the Second Piola-Kirchhoff stress tensor and the Viscosity term. */

  C = (F^(T))*F;                                 // Right Cauchy-Green strain tensor     : unitless
  J = Determinant(F);                            // J is det of F                        : unitless

  S = mu0*I + (-mu0 + 2*Lame*log(J))*(C^(-1));                                 //        : Mpa

  /* Calculate viscosity tensor:
  To do this, we need to calculate the deformation gradient. Luckily, at this
  point, the F member of P_In (the deformation gradient that P_In contains) has
  not been updated. Thus, P_In's F member is the old strain tensor. Therefore,
  we can use the F that we calculated above as the 'new' deformation tensor,
  F(t), and P_In.F as the 'old' deformation tensor, F(t-dt). We can then use
  the forward difference approximation of the derivative to get an approximation
  for F_Prine. */
  F_Prime = (1/dt)*(F - P_In.F);                                               //        : s^-1
  L = F_Prime*(F^(-1));                                                        //        : s^-1
  Visc = J*mu*(L + (L^(T))*(F^(-T)));                                          //        : Mpa

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
  Vector rj;                                     // Displacement vector                  : mm

  /* P_In aliases (notice, P_In/P_i doesn't change so these variables are const) */
  const unsigned int Num_Neighbors = P_In.Num_Neighbors;  // Number of neighbors of P_In
  const double Vi = P_In.Vol;                    // Volume of P_In                       : mm^3
  const Tensor P_i = P_In.P;                     // First Piola-Kirchhoff stress tensor  : Mpa
  const Vector * R = P_In.R;                     // Reference displacement array         : mm
  const double * W = P_In.W;                     // Shape function array                 : unitless
  const Vector * Grad_W = P_In.Grad_W;           // Grad_W array                         : 1/mm
  const double Mass = P_In.Mass;                 // P_i's mass                           : g
  const double alpha = P_In.alpha;               // alpha static member                  : unitless
  const double E = P_In.E;                       // Hourglass stiffness                  : Mpa

  /* Hour glass variables */
  double Mag_Rj;                                                               //        : mm
  double Mag_rj;                                                               //        : mm
  Vector Error_ij;                                                             //        : mm
  Vector Error_ji;                                                             //        : mm
  double delta_ij;                                                             //        : mm
  double delta_ji;                                                             //        : mm

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    // Update Neighbor
    Neighbor_ID = P_In.Neighbor_IDs[j];

    /* Calculate Internal force */
    Vj = Particles[Neighbor_ID].Vol;                                           //        : mm^3
    P_j = Particles[Neighbor_ID].P;                                            //        : Mpa
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm

    // Force_Int += Vi*Vj*(P_i + P_j)*P_In.Grad_W_Tilde[j];
    Force_Int += (Vi*Vj)*(P_i + P_j)*Grad_W[j];                                  //        : N

    /* Calculate external Force */

    /* Calculate Hour Glass force */
    Mag_Rj = Magnitude(R[j]);                                                  //        : mm
    Mag_rj = Magnitude(rj);                                                    //        : mm

    // Error_ij = F_i*R_ij - r_ij
    Error_ij = (P_In.F)*(R[j]) - rj;                                           //        : mm
    /* Note, to calculate Error_ji, we need to know the the jth particle's
    deformation gradient, R_ji, and r_ji. We could calculate all of these, but we
    can save some time by making a few clever observations. From the definion of
    R,
             R_ji = X_i - X_j = -(X_j - X_i) = -R_ij
     Likewise, from the defintion of r, we can decduce that r_ij = -r_ji. Thus,
            Error_ji = F_j*R_ji - r_ji = F_j*(-R_ij) + r_ij = -F_j*(R_ij) + r_ij
        */
    Error_ji = -1*(Particles[Neighbor_ID].F)*(R[j]) + rj;                      //        : mm

    delta_ij = Vector_Dot_Product(Error_ij, rj)/(Mag_rj);                      //        : mm
    /* Note: To find delta ji, we wan to take the dot product of Error_ji with
    r_ji. However, there's no good way of finding r_ji since that would require
    knowing the appropiate neighbor ID for particle j (Let particle i's jth
    neighbor be particle j. Particle j will have particle i as a neighbor, but
    it won't necessairally be j's 'jth' neighbor.

    Therefore, we need to get r_ji from quantities available in the ith particle.
    However, r_ji = r_i - r_j = -(r_j - r_i) = -r_ij = -rj. Thus, we can just
    use the negative of P_In.rj for r_ji! */
    delta_ji = -1*Vector_Dot_Product(Error_ji, rj)/(Mag_rj);                   //        : mm

    Force_Hg += (-1*alpha*((E*Vi*Vj*W[j])/(2*Mag_Rj*Mag_Rj*Mag_rj))*
                (delta_ij + delta_ji))*(rj);                                   //        : N
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

  /* Compute acceleration of particle at new position a(t_i+1).
  Note that all the forces we have calculated have been in units of Newtons.
  Our mass is in units of grams and we want the acceleration in units of
  mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
  multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
  then we get acceleration in mm/s^2. */
  acceleration = (1e+6)*(1./Mass)*(Force_Int + Force_Ext + Force_Hg);          //         : mm/s^2

  /* Now update the velocity, position vectors. This is done using the
  'leap-frog' integration scheme. However, during the first step of this
  scheme, we need to use forward euler to get the initial velocity.*/
  if(P_In.First_Iteration == true) {
    P_In.First_Iteration = false;
    P_In.vel = P_In.vel + (dt/2.)*acceleration;  // velocity starts at t_i+1/2           : mm/s
  } //   if(P_In.First_Iteration == true) {

  P_In.x = P_In.x + dt*P_In.vel;                 // x_i+1 = x_i + dt*v_(i+1/2)           : mm
  P_In.vel = P_In.vel + dt*acceleration;         // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s
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

#endif
