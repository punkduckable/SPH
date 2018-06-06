#if !defined(_PARTICLE_SOURCE)
#define _PARTICLE_SOURCE

////////////////////////////////////////////////////////////////////////////////
// Constructors

Particle::Particle(void) {
} // Particle::Particle(void) {

Particle::Particle(Vector X_In) {
  X = X_In;     // Set reference position
  x = X_In;     // Set current position
} // Particle::Particle(Vector X_In) {



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

void Particle::Set_Neighbor_List(const unsigned int *List,const unsigned int N) {
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

  // Calculate Grad_W using Grad_W(Rj,h) = 3*(h - |Rj|)^2*(Rj / |Rj|)
  Grad_W = 3*((h - Mag_Rj)*(h - Mag_Rj)/(Mag_Rj))*Rj;

  return Grad_W;
} // Vector Particle::Grad_W(const Vector Rj, cosnt double h) const {

#endif
