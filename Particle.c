#if !defined(_PARTICLE_SOURCE)
#define _PARTICLE_SOURCE

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  Neighbors_Set = false;
  Num_Neighbors = 0;
  V = 0;
} // Particle::Particle(void) {

Particle::Particle(const Particle & P_In) {
  /* Copy elements from P_In to local elements. Note that an object of one class
  is able to access private data of another object of that class.

  Note: because the Particle class contains pointers, we need to perform a deep
  copy of those pointers to ensure that the pointed too location isn't deleted
  when a temporary object is created and then destroyed by the copy constructor.
  */

  // Element wise copy of NON-POINTER members
  V = P_In.V;
  X = P_In.X;
  x = P_In.x;
  x_new = P_In.x_new;
  P = P_In.P;
  Neighbors_Set = P_In.Neighbors_Set;

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
  Num_Neighbors = P_In.Num_Neighbors;

  Neighbor_List = new unsigned int[Num_Neighbors];
  r = new Vector[Num_Neighbors];
  R = new Vector[Num_Neighbors];
  Grad_W = new Vector[Num_Neighbors];
  Grad_W_Tilde = new Vector[Num_Neighbors];

  for(int j = 0; j < Num_Neighbors; j++) {
    Neighbor_List[j] = P_In.Neighbor_List[j];
    r[j] = P_In.r[j];
    R[j] = P_In.R[j];
    Grad_W[j] = P_In.Grad_W[j];
    Grad_W_Tilde[j] = P_In.Grad_W_Tilde[j];
  } // for(int j = 0; j < Num_Neighbors; j++) {
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  delete [] r;
  delete [] R;
  delete [] Grad_W;
  delete [] Grad_W_Tilde;
  delete [] Neighbor_List;
} // Particle::~Particle(void) {



////////////////////////////////////////////////////////////////////////////////
// Particle setup methods:
// Set Neighbors, Set_V, Set_X.

void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_Id_List, const Particle *Particles) {
  /* First check if this particle already has neighbors. This function should
  only be called if the neighbors have not been set. The reason for this is
  that this method allocates pointers. If the pointers have already been set,
  then allocating them again will cause a memory leak.  */
  if(Neighbors_Set == true) {
    printf("Neigbors already set, returning.\n");
    return;
  }

  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Neighbor_List, and Grad_W_Tilde array
  Neighbor_List = new unsigned int[N];
  r = new Vector[N];
  R = new Vector[N];
  Grad_W = new Vector[N];
  Grad_W_Tilde = new Vector[N];

  /* Now that we know our neighbors (and where they are), there's a lot that
  we can find out. We an set the Neighbor_List, r, R, Grad_W members and
  determine the Shape tensor A (and its inverse) in order to populate the
  Grad_W_Tilde array! */

  int Neighbor_Id;                                 // Keep track of current particle
  double Mag_Rj;                                   // magnitude of Rj vector
  Tensor A, A_Inv;                                 // Shape Tensor and its inverse
  Particle P_Neighbor;                             // Current neighbor particle

  for(int j = 0; j < N; j++) {
    Neighbor_Id = Neighbor_Id_List[j];             // Get Neighbor ID (index in Particles array)
    Neighbor_List[j] = Neighbor_Id;                // Set Neighbor_List member Element
    P_Neighbor = Particles[Neighbor_Id];           // Set P_Particle to current particle (using particle's array)

    r[j] = P_Neighbor.x - x;                       // Determine spacial displacement vector
    R[j] = P_Neighbor.X - X;                       // Determine reference displacement vector

    Mag_Rj = R[j].Magnitude();                     // Determine magniude of displacement vector with jth particle (R[j])

    Grad_W[j] = 3*((h - Mag_Rj)*(h - Mag_Rj)/(Mag_Rj))*R[j];  // calculate Grad_W at jth particle

    A += Dyadic_Product(P_Neighbor.V*Grad_W[j], R[j]);   // Add in the Current Neighbor's contribution to the Shape tensor
  } // for(int j = 0; j < N; j++) {

  // Now we can calculate A^(-1) from A.
  A_Inv = A.Inverse();

  // Now we can popuate the Grad_W_Tilde array
  for(int j = 0; j < N; j++) {
    Grad_W_Tilde[j] = A_Inv*Grad_W[j];
  } // for(int j = 0; j < N; j++) {

  // Now that neighbors have been set, we set 'Neighbors_Set' to true
  Neighbors_Set = true;
} // void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_Id_List, const Particle *Particles) {

void Particle::Set_X(const Vector & X_In) {
  X = X_In;
} // void Particle::Set_X(cosnt Vector X_In) {

void Particle::Set_V(const double V_In) {
  V = V_In;
} // void Set_V(const double V_In) {

void Particle::Set_P(const Tensor & P_In) {
  P = P_In;
} // void Particle::Set_P(const Tensor & P_In) {

void Particle::Set_F(const Tensor & F_In) {
  F = F_In;
} // void Particle::Set_P(const Tensor & P_In) {


////////////////////////////////////////////////////////////////////////////////
// Operator overloading

Particle & Particle::operator=(const Particle & P_In) {
  /* Copy elements from P_In to local elements. Note that an object of one class
  is able to access private data of another object of that class.

  Note: because the Particle class contains pointers, we need to perform a deep
  copy of those pointers to ensure that the pointed too location isn't deleted
  when a temporary object is created and then destroyed by the copy constructor.
  */

  // Element wise copy of NON-POINTER members
  V = P_In.V;
  X = P_In.X;
  x = P_In.x;
  x_new = P_In.x_new;
  P = P_In.P;
  Neighbors_Set = P_In.Neighbors_Set;

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
  Num_Neighbors = P_In.Num_Neighbors;
  Neighbor_List = new unsigned int[Num_Neighbors];
  Grad_W = new Vector[Num_Neighbors];
  Grad_W_Tilde = new Vector[Num_Neighbors];

  for(int i = 0; i < Num_Neighbors; i++) {
    Neighbor_List[i] = P_In.Neighbor_List[i];
    Grad_W[i] = P_In.Grad_W[i];
    Grad_W_Tilde[i] = P_In.Grad_W_Tilde[i];
  } // for(int i = 0; i < Num_Neighbors; i++) {

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
// Calculate W

double Particle::Calc_W(const Vector & Rj) {
  double Mag_Rj = Magnitude(Rj);

  // If vector is outside of support radius, return 0
  if(Mag_Rj > h)
    return 0;

  return (h - Mag_Rj)*(h - Mag_Rj)*(h - Mag_Rj);
} // double Calc_W(const Vector & Rj_In) {

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
  these arguments to calculate P (the first Piola-Kirchoff stress tensor) */

  /* First, let's set up the local variables that will be used to update the
  particle's position */
  double Vj;                                          // Volume of jth particle (a neighbor)
  int Neighbor_Id;                                    // Index of jth neighbor.

  Tensor F_New;                                       // New Deformation gradient

  Tensor C;                                           // Richt-Cauchy stress tensor
  Tensor S;                                           // Second Poila-Kirchoff stress tensor
  Vector M;                                           // Fber orientation vector
  Tensor M_Dyad_M;                                    // Stores dyadic product of M with M... used to calculate strain energy
  Tensor I = {1,0,0,
              0,1,0,
              0,0,1};                                 // Identity tensor

  double k1, k2, mu0, mu;
  double I1, J, I4, p;                                // Used to calculate stress energy function

  Tensor F_Prime;                                     // Approximates the time derivative of F
  Tensor L;                                           // Used in viscosity correction term
  Tensor Visc;                                        // Viscosity correction term for P


  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) Vj Grad_W(Rj, h) */
  for(int j = 0; j < P_In.Num_Neighbors; j++) {
    Neighbor_Id = P_In.Neighbor_List[j];
    Vj = Particles[Neighbor_Id].V;

    F_New += Dyadic_Product(P_In.r[j], Vj*P_In.Grad_W_Tilde[j]);
  } // for(int j = 0; j < Num_Neighbors; j++) {

  /* Now that we have calculated the deformation gradient, we need to calculate
  the first Piola-Kirchhoff stess tensor. To do this, however, we need to
  find the Second Piola-Kirchhoff stress tensor and the Viscosity term. */

  // Calculate S
  k1 = P_In.k1;
  k2 = P_In.k2;
  mu0 = P_In.mu0;
  mu = P_In.mu;
  M = P_In.M;

  C = F_New.Transpose()*F_New;
  I1 = C(1,1) + C(2,2) + C(3,3);
  J = Determinant(F_New);
  M_Dyad_M = Dyadic_Product(M,M);
  I4 = Tensor_Dot_Product( C, M_Dyad_M );
  p = mu0 / (C(1,1)*C(2,2) - C(1,2)*C(2,1));

  S = mu0*I + k1*(I4 - 1)*exp(k1*(I4-1)*(I4-1))*M_Dyad_M - p*(1./J)*Inverse(C);

  // Calculate viscosity tensor
  F_Prime = (1/dt)*(F_New - P_In.F);
  L = F_Prime*Inverse(F_New);
  Visc = J*mu*(L + L.Transpose())*Transpose(F_New.Inverse());

  // Set P
  P_In.Set_P(F_New*S + Visc);

  // Update Particle's F member
  P_In.Set_F(F_New);

} // void Update_P(const Particle & P_In, const Particle * Particles, const double dt) {

void Update_Particle_Position(Particle & P_In, const Particle * Particles, const double dt) {
  /* This function assumes that every particle in the Particle's array has
  an updated P tensor. Likewise, it assumes that the E and alpha static
  member variables have been set. This function should not be run until
  these assumptions are valid. */

  Vector Force_Int;                                  // Internal Force vector
  Vector Force_Ext;                                  // External/body force
  Vector Force_Hg;                                   // Hour-glass force

  Vector acceleration;                               // acceleration vector

  int Neighbor_Id;                                   // ID of current neighbor particle (in paritlce's array)
  int Num_Neighbors = P_In.Num_Neighbors;            // Number of neighbors of P_In.
  Particle P_Neighbor;                               // Current neighbor particle


  double Vj;                                         // Volume of jth particle
  double Vi = P_In.V;                                // Volume of P_In;
  Tensor P_j;                                         // Piola-Kirchhoff stress tensor for jth particle
  Tensor P_i = P_In.P;                                // Piola-Kirchhoff stress tensor for P_In

  Vector Rj;
  Vector rj;
  double Mag_Rj;
  double Mag_rj;
  Vector Error_ij;
  Vector Error_ji;
  double delta_ij;
  double delta_ji;

  double alpha = P_In.alpha;
  double E = P_In.E;
  double rho = P_In.rho;


  for(int j = 0; j < Num_Neighbors; j++) {
    // Update Neighbor
    Neighbor_Id = P_In.Neighbor_List[j];
    P_Neighbor = Particles[Neighbor_Id];

    // Calculate Internal force
    Vj = P_Neighbor.V;

    P_j = P_Neighbor.P;

    Force_Int += Vi*Vj*(P_i + P_j)*P_In.Grad_W_Tilde[j];

    // Calculate external Force

    // Calculate Hour Glass force
    Rj = P_In.R[j];
    Mag_Rj = Magnitude(Rj);
    rj = P_In.r[j];
    Mag_rj = Magnitude(rj);

    Error_ij = (P_In.F)*(Rj) - rj;
    Error_ji = (P_Neighbor.F)*(-1*Rj) + rj;                    // Rij = -Rji and rij = -rij

    delta_ij = Vector_Dot_Product(Error_ij, P_In.r[j])/(Mag_rj);
    delta_ji = Vector_Dot_Product(Error_ji, P_Neighbor.r[j])/(Mag_rj);

    Force_Hg += -1*((alpha*E*Vi*Vj*P_In.Calc_W(Rj))/(2*Mag_Rj*Mag_Rj))*(delta_ij + delta_ji)*(rj/Mag_rj);
  } // for(int j = 0; j < Num_Neighbors; j++) {

  // Compute acceleration
  acceleration = (1./(Vi*rho))*(Force_Int + Force_Ext + Force_Hg);

  // Calculate new position
}

#endif
