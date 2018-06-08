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
// Set neighbor List

void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_List_In, const Vector *X_Neighbors) {
  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Neighbor_List, and Grad_W_Tilde array
  Neighbor_List = new int[N];
  Grad_W_Tilde = new Vector[N];

  // Calculate R from X_Neibhors. Rj = Xj - X (where is is current particle, j is some neighbor particle)
  Vector R[N];
  for(int j = 0; j < N; j++) {
    R[j] = X_Neighbors[j] - X;
  }

  // Now that we have our list of neighbors, we can populate Populate neighbor list using supplied list
  for(int j = 0; j < N; j++) {
    Neighbor_List[j] = Neighbor_List_In[j];
    Grad_W_Tilde[j] = Calc_Grad_W_Tilde(R[j]);
  } // for(int j = 0; j < N; j++) {
} //void Particle::Set_Neighbor_List(const unsigned int *List,const unsigned int N) {



////////////////////////////////////////////////////////////////////////////////
// Update x (spacial position)


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

void Particle::Calc_A_Inv(void) {
  /* this function is used to calculate the inverse of the shape tensor A.
  To do this, I need to first find the shape tensor A. This function should only
  be run once as soon as the particle is setup. The shape matrix is defined by,
  A = sum of all j in S { Vj*Grad_W(rj,h) (dyadic product) Rj }
  Once we have found this tensor, we can find its inverse using the inverse
  method for tensors.

  Note, when we initialize the tensor A, the default tensor constructor sets
  every element to 0. */

  int Neighbor_Id;
  Vector Rj;
  Tensor A;

  for(int i = 0; i < Num_Neighbors; i++) {
    Neighbor_Id = Neighbor_List[i];
    Rj = Particles[Neighbor_Id].Get_X();

    A += Dyadic_Product(Vj*( Grad_W(Rj)) , Rj );
  } // for(int i = 0; i < Num_Neighbors; i++) {

  A_Inv = A.Inverse();
} // void Particle::Calculate_A_Inv() {

Vector Particle::Calc_Grad_W_Tilde(const Vector Rj) const {
  return A_Inv*(Grad_W(Rj));
} // Vector Particle::Calculate_Grad_W_Tilde(const Vector Rj) const {

#endif
