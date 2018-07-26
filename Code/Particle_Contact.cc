#if !defined(PARTICLE_CONTACT)
#define PARTICLE_CONTACT

#include "Particle_Helpers.h"

void Particle_Helpers::Contact(Particle_Array & Body_A, Particle_Array & Body_B) {
  // Get paramaters from Body_A, Body_B
  const unsigned int Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned int Num_Particles_B = Body_B.Get_Num_Particles();
  const double Shape_Function_Amp = Body_A.Get_Shape_Function_Amplitude();

  /* This function implements Particle-Particle 'contact'.
  the contact force is applied whenever two particles are within a 2h radius of
  one another (when the two particles's support radii overlap). The applied
  force is in the direction of the line between the two particles centers. */
  const double h = (Body_A.Get_h() + Body_B.Get_h())/2;    // Contact force support radius         : mm
  const double h_squared = h*h;                            // Square of support radius   : mm^2
  unsigned int i,j;
  double V_i, V_j;                                                             //        : mm
  double K_V_i;                                                                //        : N*mm^2
  double Mag_r_ij;                                                             //        : mm
  Vector r_ij, Grad_W;                                                         //        : mm Vector
  Vector * Body_B_x = new Vector[Num_Particles_B];                             //        : mm Vector
  Vector F_Contact;                                                            //        : N Vector

  // First, set every particle's Contact force to zero. Also, store all of
  // Body B's position vectors in an array (reduces the number of cache misses)
  for(i = 0; i < Num_Particles_A; i++)
    Body_A[i].Force_Contact = {0,0,0};                                         //        : N Vector

  for(i = 0; i < Num_Particles_B; i++) {
    Body_B[i].Force_Contact = {0,0,0};                                         //        : N Vector
    Body_B_x[i] = Body_B[i].x;                                                 //        : mm Vector
  }

  // For each particle in A, check if there is a particle in B that is within
  // the combined support radius
  for(i = 0; i < Num_Particles_A; i++) {
    V_i = Body_A[i].Vol;                                                       //        : mm^3
    K_V_i = K*V_i;                                                             //        : N*mm

    for(j = 0; j < Num_Particles_B; j++) {
      V_j = Body_B[j].Vol;                                                     //        : mm^3
      r_ij = Body_A[i].x - Body_B_x[j];                                        //        : mm Vector

      // Check if |rij| < h. Note that this is equivalent to rij dot rij < h^2
      // If so, add contact force
      if(Vector_Dot_Product(r_ij, r_ij) < h_squared) {
        // Calculate the contact force
        Mag_r_ij = Magnitude(r_ij);                                            //        : mm
        Grad_W = (-3*(Shape_Function_Amp)*((h - Mag_r_ij)*(h - Mag_r_ij))/Mag_r_ij)*(r_ij);    // 1/(mm^4) Vector
        F_Contact = (K_V_i*V_j)*Grad_W;                                        //        : N Vector

        // Now apply the force to the two interacting bodies (Note the forces
        // are equal and opposite)
        Body_A[i].Force_Contact -= F_Contact;                                  //        : N Vector
        Body_B[j].Force_Contact += F_Contact;                                  //        : N Vector
      } // if(Magnitude(r_ij) < h) {
    } // for(j = 0; j < Num_Particle_B, j++) {
  } // for(i = 0; i < Num_Particles_A; i++) {

  delete [] Body_B_x;                                                          //        : mm Vector
} // void Particle_Helpers::Contact(Particle_Array & Body_A, Particle_Array & Body_B) {

#endif
