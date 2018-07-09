#if !defined(PARTICLE_CONTACT)
#define PARTICLE_CONTACT

#include "Particle_Helpers.h"

void Particle_Helpers::Contact(Particle * Body_A, const unsigned int Num_Particles_A,Particle * Body_B, const unsigned int Num_Particles_B, const double h_In) {
  /* This function implements Particle-Particle 'contact'.
  the contact force is applied whenever two particles are within a 2h radius of
  one another (when the two particles's support radii overlap). The applied
  force is in the direction of the line between the two particles centers. */
  const double h = 2*h_In;                      // Contact force support radius
  const double h_squared = h*h;                 // Square of support radius (for performance)
  unsigned int i,j;
  double V_i, V_j, K_V_i;
  double Mag_r_ij;
  Vector r_ij, Grad_W;
  Vector * Body_B_x = new Vector[Num_Particles_B];

  // First, set every particle's Contact force to zero. Also, store all of
  // Body B's position vectors in an array (reduces the number of cache misses)
  for(i = 0; i < Num_Particles_A; i++)
    Body_A[i].Force_Contact = {0,0,0};

  for(i = 0; i < Num_Particles_B; i++) {
    Body_B[i].Force_Contact = {0,0,0};
    Body_B_x[i] = Body_B[i].x;
  }

  // For each particle in A, check if there is a particle in B that is within
  // the combined support radius
  for(i = 0; i < Num_Particles_A; i++) {
    V_i = Body_A[i].Vol;
    K_V_i = K*V_i;

    for(j = 0; j < Num_Particles_B; j++) {
      V_j = Body_B[j].Vol;
      r_ij = Body_A[i].x - Body_B_x[j];

      // Check if |rij| < h. Note that this is equivalent to rij dot rij < h^2
      // If so, add contact force
      if(Vector_Dot_Product(r_ij, r_ij) < h_squared) {
        Mag_r_ij = Magnitude(r_ij);
        Grad_W = (-3*(Particle::Shape_Function_Amp/64.)*((h - Mag_r_ij)*(h - Mag_r_ij))/Mag_r_ij)*(r_ij);
        Body_A[i].Force_Contact -= (K_V_i*V_j)*Grad_W;
        Body_B[j].Force_Contact += Body_A[i].Force_Contact;
      } // if(Magnitude(r_ij) < h) {
    } // for(j = 0; j < Num_Particle_B, j++) {
  } // for(i = 0; i < Num_Particles_A; i++) {

  delete [] Body_B_x;
} // void Particle_Helpers::Contact(Particle * Body_A, const unsigned int Num_Particles_A,Particle * Body_B, const unsigned int Num_Particles_B, const double h_In) {

#endif
