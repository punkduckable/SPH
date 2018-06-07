#if !defined(_PARTICLE_SOURCE)
#define _PARTICLE_SOURCE

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  ++Num_Particles;
} // Particle::Particle(void) {

Particle::Particle(const Vector & X_In) {
  ++Num_Particles;
  X = X_In;     // Set reference position
  x = X_In;     // Set current position
} // Particle::Particle(const Vector & X_In) {

Particle::Particle(const Particle & P_In) {
  /* Copy elements from P_In to local elements. Note that an object of one class
  is able to access private data of another object of that class. */

  ++Num_Particles;

  V = P_In.V;
  X = P_In.X;
  x = P_In.x;
  Force = P_In.Force;
  Acceleration = P_In.acceleration;
  F = P_In.F;
  A_Inv = P_In.A_Inv;
  P = P_In.P;
  S = P_In.S;
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  --Num_Particles;
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

void Particle::Set_Neighbor_List(const unsigned int *List, const unsigned int N) {
  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Neighbor_List
  Neighbor_List = (unsigned int *)malloc(N*sizeof(unsigned int));

  // Populate neighbor list using supplied list
  for(int i = 0; i < N; i++) {
    Neighbor_List[i] = List[i];
  } //   for(int i = 0; i < N; i++) {
} //void Particle::Set_Neighbor_List(const unsigned int *List,const unsigned int N) {



////////////////////////////////////////////////////////////////////////////////
// Update x (spacial position)


Vector Particle::Grad_W(const Vector Rj,const double h) const {
  // Set up a vector to hold the gradient at the jth particle.
  Vector Grad_W;

  // Calculate the magnitude of Rj
  double Mag_Rj = Rj.Magnitude();

  /* Calculate Grad_W using Grad_W(Rj,h) = (d/d|Rj|)W(Rj,h)*(Rj / |Rj|) with
  W(Rj,h) = (h - |Rj|)^3. Thus,
  Grad_W(Rj,h) = 3*(h - |Rj|)^2*(Rj / |Rj|) */
  Grad_W = 3*((h - Mag_Rj)*(h - Mag_Rj)/(Mag_Rj))*Rj;

  return Grad_W;
} // Vector Particle::Grad_W(const Vector Rj, cosnt double h) const {

void Particle::Calculate_A_Inv(const double h) {
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

    A += Dyadic_Product(Vj*( Grad_W(Rj,h)) , Rj );
  } // for(int i = 0; i < Num_Neighbors; i++) {

  A_Inv = A.Inverse();
} // void Particle::Calculate_A_Inv(const double h) {

Vector Particle::Grad_W_Tilde(const Vector Rj,const double h) const {
  return A_Inv*(Grad_W(Rj,h));
} // Vector Particle::Calculate_Grad_W_Tilde(const Vector Rj,const double h) const {

#endif
