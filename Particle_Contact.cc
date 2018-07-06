#if !defined(PARTICLE_CONTACT)
#define PARTICLE_CONTACT

#include "Particle_Helpers.h"

void Particle_Helpers::Contact(Particle * Body_A, const unsigned int Num_Particles_A,Particle * Body_B, const unsigned int Num_Particles_B, const double h_particle, const double IPS) {
  const double h = 2*h_particle;
  const double h_squared = h*h;
  unsigned int i,j;
  double V_i, V_j;
  double Mag_r_ij;
  Vector r_ij, Grad_W;

  // First, set every particle's Contact force to zero.
  for(i = 0; i < Num_Particles_A; i++)
    Body_A[i].Force_Contact = {0,0,0};

  for(i = 0; i < Num_Particles_B; i++)
    Body_B[i].Force_Contact = {0,0,0};

  // For each particle in A, check if there is a particle in B that is within
  // the combined support radius
  for(i = 0; i < Num_Particles_A; i++) {
    V_i = Body_A[i].Vol;

    for(j = 0; j < Num_Particles_B; j++) {
      V_j = Body_B[j].Vol;
      r_ij = Body_A[i].x - Body_B[j].x;

      // Check if |rij| < h. If so, implement a contact force
      if(Vector_Dot_Product(r_ij, r_ij) < h_squared) {
        Mag_r_ij = Magnitude(r_ij);
        Grad_W = -3*(Particle::Shape_Function_Amp/64.)*((h - Mag_r_ij)*(h - Mag_r_ij))*(r_ij / Mag_r_ij);
        Body_A[i].Force_Contact -= V_i*V_j*K*Grad_W;
        Body_B[j].Force_Contact += Body_A[i].Force_Contact;
      } // if(Magnitude(r_ij) < h) {
    } // for(j = 0; j < Num_Particle_B, j++) {
  } // for(i = 0; i < Num_Particles_A; i++) {
} // void Particle_Helpers::Contact(Particle * Body_A, const unsigned int Num_Particles_A,Particle * Body_B, const unsigned int Num_Particles_B, const double h_particle, const double IPS) {

#endif
