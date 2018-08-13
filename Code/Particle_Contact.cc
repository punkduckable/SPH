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
  const double h = (Body_A.Get_Inter_Particle_Spacing() + Body_B.Get_Inter_Particle_Spacing());    // Contact force support radius         : mm
  const double h_squared = h*h;                            // Square of support radius   : mm^2
  unsigned int i,j;
  double V_i, V_j;                                                             //        : mm
  double K_V_i;                                                                //        : N*mm^2
  double Mag_r_ij;                                                             //        : mm
  Vector r_ij, Grad_W;                                                         //        : mm Vector
  Vector * Body_B_x = new Vector[Num_Particles_B];                             //        : mm Vector
  Vector F_Contact;                                                            //        : N Vector

  // Thread-local (private) force contributions to Body_B (see description below)
  Vector * Body_B_F_Contact_Local = new Vector[Num_Particles_B];               //        : N Vector

  // First, store all of Body B's position vectors in an array (this improves performance... I think?)
  for(i = 0; i < Num_Particles_B; i++) {
    Body_B_x[i] = Body_B[i].Get_x();                                           //        : mm Vector
    Body_B_F_Contact_Local[i] = {0,0,0};
  }

  // For each particle in A, check if there is a particle in B that is within
  // the combined support radius
  #pragma omp for schedule(dynamic)
  for(i = 0; i < Num_Particles_A; i++) {
    // Skip broken particles
    if(Body_A[i].Get_D() >= 1)
      continue;

    V_i = Body_A[i].Get_Vol();                                                 //        : mm^3
    K_V_i = K*V_i;                                                             //        : N*mm

    for(j = 0; j < Num_Particles_B; j++) {
      V_j = Body_B[j].Get_Vol();                                               //        : mm^3
      r_ij = Body_A[i].Get_x() - Body_B_x[j];                                  //        : mm Vector

      // Check if |rij| < h. Note that this is equivalent to rij dot rij < h^2
      // If so, add contact force
      if(Vector_Dot_Product(r_ij, r_ij) < h_squared) {
        // Calculate the contact force
        Mag_r_ij = Magnitude(r_ij);                                            //        : mm
        Grad_W = (-3*(Shape_Function_Amp)*((h - Mag_r_ij)*(h - Mag_r_ij))/Mag_r_ij)*(r_ij);    // 1/(mm^4) Vector
        F_Contact = (K_V_i*V_j)*Grad_W;                                        //        : N Vector

        /* Now apply the force to the two interacting bodies (Note the forces
        are equal and opposite). It should be noted that we don't actually
        apply the force to body B's jth particle. The reason for this is
        race conditions. It is possible for two threads to attempt to write
        the contact force to the same particle in body B at the same time.
        This causes a race condition (when run in parallel). This can be
        fixed using a critical region, but that's slow. To fix this, each
        thread gets its own 'Local contact force array'. each thread stores
        the force that it would apply to the particles in body B in this array.
        Once we have cycled through all of body A's particles, we then 'reduce'
        these force arrays. One by one, the threads add their contribuitions to
        body_b's contact forces. This runs about 2x as quick as if we had
        inserted the critical region into this loop. */
        Body_A[i].Force_Contact -= F_Contact;                                  //        : N Vector
        Body_B_F_Contact_Local[j] += F_Contact;                                //        : N Vector
      } // if(Magnitude(r_ij) < h) {
    } // for(j = 0; j < Num_Particle_B, j++) {
  } // for(i = 0; i < Num_Particles_A; i++) {

  /* Now that we have finished the loop above, each thread has its
  contribution to the contact force for each of Body_B's particles. We therefore
  add these contributions to Body_B's particles one by one (using a critical
  region) */
  #pragma omp critical
  for(j = 0; j < Num_Particles_B; j++)
    Body_B[j].Force_Contact += Body_B_F_Contact_Local[j];                      //        : N Vector

  delete [] Body_B_x;                                                          //        : mm Vector
} // void Particle_Helpers::Contact(Particle_Array & Body_A, Particle_Array & Body_B) {

#endif
