#include "Body.h"

void Body::Contact(Body & Body_A, Body & Body_B) {
  // Get paramaters from Body_A, Body_B
  const unsigned Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned Num_Particles_B = Body_B.Get_Num_Particles();
  const double Shape_Function_Amp = Body_A.Get_Shape_Function_Amplitude();

  /* This function implements Particle-Particle 'contact'.
  the contact force is applied whenever two particles are within a 2h radius of
  one another (when the two particles's support radii overlap). The applied
  force is in the direction of the line between the two particles centers. */
  const double h = Simulation::Contact_Distance;           // Contact force support radius         : mm
  const double h_squared = h*h;                            // Square of support radius   : mm^2
  unsigned i;
  Vector r_ij;                                                                 //        : mm Vector
  Vector * Body_B_x = new Vector[Num_Particles_B];                             //        : mm Vector
  Vector Grad_W;                                                               //        : 1/mm^4 Vector
  Vector x_i;                                                                  //        : mm Vector
  Vector F_Contact, F_Friction;                                                //        : N Vector
  Vector Relative_Velocity;                                                    //        : mm/s Vector
  bool Contact_Flag = 0;

  // Thread-local (private) force contributions to Body_B (see description below)
  Vector * Body_B_F_Contact_Local = new Vector[Num_Particles_B];               //        : N Vector
  Vector * Body_B_F_Friction_Local = new Vector[Num_Particles_B];              //        : N Vector

  // First, store all of Body B's position vectors in an array (this improves performance... I think?)
  for(i = 0; i < Num_Particles_B; i++) {
    Body_B_x[i] = Body_B[i].Get_x();                                           //        : mm Vector
    Body_B_F_Contact_Local[i] = {0,0,0};                                       //        : N Vector
    Body_B_F_Friction_Local[i] = {0,0,0};                                      //        : N Vector
  } //   for(i = 0; i < Num_Particles_B; i++) {

  // For each particle in A, check if there is a particle in B that is within
  // the combined support radius
  #pragma omp for schedule(dynamic)
  for(i = 0; i < Num_Particles_A; i++) {
    // Skip broken particles
    if(Body_A[i].Get_D() >= 1) { continue; }

    double V_i = Body_A[i].Get_Vol();                                          //        : mm^3
    const double KV_i = K*V_i;                                                 //        : N*mm
    x_i = Body_A[i].Get_x();                                                   //        : mm Vector

    for(unsigned j = 0; j < Num_Particles_B; j++) {
      r_ij = x_i - Body_B_x[j];                                                //        : mm Vector

      // Check if |rij| < h. Note that this is equivalent to rij dot rij < h^2
      // If so, add contact force
      if(Dot_Product(r_ij, r_ij) < h_squared) {
        /* First, set the 'contact flag' to true. This flag is designed to
        improve perfomance. The idea here is that there will be no contact
        for most of the time steps. Because of this, it doesn't make sense
        to run throug the critical for loop (to add each thread's contact and
        friction forces to the particles). Thus, we only run that critical loop
        if this flag is false. The idea here is that the cost of updating this
        flag for each contacting particle is cheaper than running through the
        critical for loop on every time step. */
        Contact_Flag = true;

        // Calculate the contact force
        double V_j = Body_B[j].Get_Vol();                                      //        : mm^3
        double Mag_r_ij = Magnitude(r_ij);                                     //        : mm
        double h_minus_Mag_r_ij = h - Mag_r_ij;                                //        : mm
        Grad_W = (-3*(Shape_Function_Amp)*(h_minus_Mag_r_ij*h_minus_Mag_r_ij)/Mag_r_ij)*(r_ij);   // 1/mm^4 Vector

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
        body_b's contact forces. This runs about 2x faster than if we had
        inserted the critical region into this loop. */
        F_Contact = (KV_i*V_j)*Grad_W;                                         //        : N Vector
        Body_A[i].Force_Contact -= F_Contact;                                  //        : N Vector
        Body_B_F_Contact_Local[j] += F_Contact;                                //        : N Vector

        /* Now let's calculate the frictional force. To do this, we first need
        to get the unit vector of the relative velocity between bodies A and B
        as well as the magnitude of the contact force */
        double Mag_F_Contact = F_Contact.Magnitude();                          //        : N
        Relative_Velocity = Body_A[i].V - Body_B[j].V;                         //        : mm/s Vector

        // Now we can calculate the frictional force and apply it to the two bodies
        F_Friction = ((-1*FRICTION_COEFFICIENT*Mag_F_Contact) / Relative_Velocity.Magnitude())*(Relative_Velocity);
        Body_A[i].Force_Friction += F_Friction;                                //        : N Vector
        Body_B_F_Friction_Local[j] -= F_Friction;                              //        : N Vector
      } // if(Magnitude(r_ij) < h) {
    } // for(j = 0; j < Num_Particle_B, j++) {
  } // for(i = 0; i < Num_Particles_A; i++) {

  /* Now that we have finished the loop above, each thread has its
  contribution to the contact force for each of Body_B's particles. We therefore
  add these contributions to Body_B's particles one by one (using a critical
  region) */
  #pragma omp critical
  if(Contact_Flag) {
    for(i = 0; i < Num_Particles_B; i++) {
      Body_B[i].Force_Contact += Body_B_F_Contact_Local[i];                      //        : N Vector
      Body_B[i].Force_Friction += Body_B_F_Friction_Local[i];                    //        : N Vector
    } // for(i = 0; i < Num_Particles_B; i++) {
  } // if(Contact_Flag) {

  // Now free any dynamically allocated memory.
  delete [] Body_B_x;
  delete [] Body_B_F_Contact_Local;
  delete [] Body_B_F_Friction_Local;

  /* Note, there is no explicit barrier here. This is because the next kernel
  , update_x , has a for loop before it uses the forces modified above. This
  means that the update_x method will force synchronize each thread before any
  data from here is used. Therefore, we don't need a barrier. */
} // void Body::Contact(Body & Body_A, Body & Body_B) {



void Body::Print_Net_External_Force(const unsigned l) {
  /* This function is used to find and print the net external force on a body
  This function can NOT be called by multiple threads at once (this
  function is not thread safe). */

  // First, open the file.
  FILE * File;
  if(Times_Printed_Net_External_Force == 0) {
    File = fopen("../Files/Force_Files/Net_External_Force.txt","w");
  } // if(Times_Printed_Net_External_Force == 0) {
  else {
    File = fopen("../Files/Force_Files/Net_External_Force.txt","a");
  } // else {

  // Increment the number of times that we're printed net force data.
  Times_Printed_Net_External_Force++;

  // Now add up net external force on supplied particle array and print it out
  // Note that we must do this using a single thread
  Vector Net_Contact_Force = {0,0,0};

  for(unsigned i = 0; i < Num_Particles; i++) {
    Net_Contact_Force += Particles[i].Get_Force_Friction();
    Net_Contact_Force += Particles[i].Get_Force_Contact();
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  fprintf(File,"%6u:  <%10.4f, %10.4f, %10.4f>\n", l, Net_Contact_Force(0), Net_Contact_Force(1), Net_Contact_Force(2));

  // Now close the file.
  fclose(File);
} // void Body::Print_Net_External_Force(const unsigned l) {
