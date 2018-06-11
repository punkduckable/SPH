#if !defined(_PARTICLE_SOURCE)
#define _PARTICLE_SOURCE

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
} // Particle::Particle(void) {

Particle::Particle(const Vector & X_In) {
  X = X_In;     // Set reference position
  x = X_In;     // Set current position
} // Particle::Particle(const Vector & X_In) {

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
  Force = P_In.Force;
  Acceleration = P_In.acceleration;
  F = P_In.F;
  A_Inv = P_In.A_Inv;
  P = P_In.P;
  S = P_In.S;

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
  Neighbor_List = new int[Num_Neighbors];
  Grad_W_Tilde = new Vector[Num_Neighbors];

  for(int i = 0; i < Num_Neighbors; i++) {
    Neighbor_list[i] = P_In.Neighbor_List[i];
    Grad_W_Tilde[i] = P_In.Grad_W_Tilde[i];
  } // for(int i = 0; i < Num_Neighbors; i++) {
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  delete [] Grad_W_Tilde;
  delete [] Neighbor_List;
} // Particle::~Particle(void) {



////////////////////////////////////////////////////////////////////////////////
// Set Neighbors + Helper methods  (Calc_Grad_W + Calc_Grad_W_Tilde)

void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_List_In, const Vector *X_Neighbors, const *double Neighbor_Volume) {
  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Neighbor_List, and Grad_W_Tilde array
  Neighbor_List = new int[N];
  Grad_W_Tilde = new Vector[N];

  /* Now that we know our neighbors (and where they are), we can find our shape
  Tensor and use this to calculate our Grad_W (Grad_W_Tilde). To begin, however,
  we must first calculate the displacement (in the reference configuration)
  between us and our neighbors. We store these positions in a vector called R. */
  Vector *R = new Vector[N];       // Memory must be dynamically allocated since N is variable
  for(int j = 0; j < N; j++) {
    R[j] = X_Neighbors[j] - X;
  }

  // Now that we have the displacements, we can calculate the shape tensor A!
  Tensor A;
  for(int j = 0; k < N; k++) {
    A += Dyadic_Product(Neighbor_Volume[j]*Grad_W(R[j]), R[j]);
  } // for(int j = 0; j < N; j++) {

  // Now we can calculate A^(-1) from A.
  A_Inv = A.Inverse();

  // Now we can popuate the Neighbor_List and Grad_W_Tilde members!
  for(int j = 0; j < N; j++) {
    Neighbor_List[j] = Neighbor_List_In[j];
    Grad_W_Tilde[j] = Calc_Grad_W_Tilde(R[j]);
  } // for(int j = 0; j < N; j++) {
} //void Particle::Set_Neighbor_List(const unsigned int *List,const unsigned int N) {

Vector Particle::Calc_Grad_W_Tilde(const Vector Rj) const {
  return A_Inv*(Grad_W(Rj));
} // Vector Particle::Calculate_Grad_W_Tilde(const Vector Rj) const {

Vector Particle::Calc_Grad_W(const Vector Rj) const {
  // Set up a vector to hold the gradient at the jth particle.
  Vector Grad_W;

  // Calculate the magnitude of Rj
  double Mag_Rj = Rj.Magnitude();

  /* Calculate Grad_W using Grad_W(Rj,h) = (d/d|Rj|)W(Rj,h)*(Rj / |Rj|) with
  W(Rj,h) = (h - |Rj|)^3. Thus,
  Grad_W(Rj,h) = 3*(h - |Rj|)^2*(Rj / |Rj|) */
  Grad_W = 3*((h - Mag_Rj)*(h - Mag_Rj)/(Mag_Rj))*Rj;

  return Grad_W;
} // Vector Particle::Grad_W(const Vector Rj) const {



////////////////////////////////////////////////////////////////////////////////
// Position functions

void Particle::Set_X(const Vector X_In) {
  X = X_In;
} // void Particle::Set_X(cosnt Vector X_In) {

Vector Particle::Get_X(void) const {
  return X;
} // Vector Particle::Get_X(void) const {

Vector Particle::Get_x(void) cosnt {
  return x;
} // Vector Particle::Get_x(void) cosnt {



////////////////////////////////////////////////////////////////////////////////
// Operator overloading

Particle & Particle::operator=(const Particle P_In) {
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
  Force = P_In.Force;
  Acceleration = P_In.acceleration;
  F = P_In.F;
  A_Inv = P_In.A_Inv;
  P = P_In.P;
  S = P_In.S;

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
  Neighbor_List = new int[Num_Neighbors];
  Grad_W_Tilde = new Vector[Num_Neighbors];

  for(int i = 0; i < Num_Neighbors; i++) {
    Neighbor_list[i] = P_In.Neighbor_List[i];
    Grad_W_Tilde[i] = P_In.Grad_W_Tilde[i];
  } // for(int i = 0; i < Num_Neighbors; i++) {

  return *this;
}



////////////////////////////////////////////////////////////////////////////////
// Friend functions (Update x)

void Update_Particle_Position(const Particle & P_In, const Particle * Particles, const double dt) {
  /* The purpose of this function is to update the spacial position of particle
  P_In. This function assumes that the position of each of P_In's neighbors has
  been updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its Grad_W_Tilde
  vector. This function should not be called until these assumptions have been
  met.

  What does this function do? It accepts a Particle (P_In), a list of all
  particles in the current body, and the desired time step. This function first
  calculates the acceleration of P_In given the current configuration, and then
  uses this to find P_In's new position. This new position is stored in the
  x_new member variable (which, using a different function, is moved to x after
  we have found x_new for every particle). */

  /* First, let's set up the local variables that will be used to update the
  particle's position */
  double Vj;                                          // Volume of jth particle (a neighbor)
  int Neighbor_Id;                                    // Index of jth neighbor.

  Vector Force;                                       // Force vector
  Vector acceleration;                                // acceleration vector
  Vector rj;                                          // Spacial displacement to jth particle (a neighbor)

  Tensor F;                                           // Deformation gradient
  Tensor A_Inv;                                       // Inverse of shape tensor
  Tensor P;                                           // First Poila-Kirchoff stress tensor
  Tensor S;                                           // Second Poila-Kirchoff stress tensor

  Particle P_Neighbor;                                // Current neighbor Particle

  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) Vj Grad_W(Rj, h) */
  for(int j = 0; j < Num_Neighbors; j++) {
    Neighbor_Id = P_In.Neighbor_List[j];
    P_Neighbor = Particles[Neighbor_Id];
    rj = P_Neighbor.x - P_In.x;

    F += Dyadic_Product(rj, P_In.Grad_W_Tilde[j]);
  }

  /* Now we can calculate P */
} // void Update_Particle_Position(const Particle & P_In, const Particle * Particles, const double dt) {

#endif
