#include "Body.h"
#include "Simulation/Simulation.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Array.h"
#include <math.h>

void Body::Contact(Body & Body_A, Body & Body_B) {
  // Get paramaters from Body_A, Body_B
  const unsigned Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned Num_Particles_B = Body_B.Get_Num_Particles();
  const double Shape_Function_Amp = Body_A.Get_Shape_Function_Amplitude();

  /* This function implements Particle-Particle 'contact'.
  the contact force is applied whenever two particles are within a 2h radius of
  one another (when the two particles's support radii overlap). The applied
  force is in the direction of the line between the two particles' centers. */
  const double h = Simulation::Contact_Distance;           // Contact force support radius         : mm
  const double h_squared = h*h;                            // Square of support radius   : mm^2
  unsigned i;
  Vector r_ij;                                                                 //        : mm Vector
  Vector * Body_B_x = new Vector[Num_Particles_B];                             //        : mm Vector
  Vector Grad_W;                                                               //        : 1/mm^4 Vector
  Vector x_i;                                                                  //        : mm Vector
  Vector F_Contact, F_Friction;                                                //        : N Vector
  Vector Relative_Velocity;                                                    //        : mm/s Vector
  bool Contact_Flag = false;

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

    double V_i = Body_A[i].Get_Volume();                                       //        : mm^3
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
        double V_j = Body_B[j].Get_Volume();                                   //        : mm^3
        double Mag_r_ij = Magnitude(r_ij);                                     //        : mm
        double h_minus_Mag_r_ij = h - Mag_r_ij;                                //        : mm
        Grad_W = (-3*(Shape_Function_Amp)*(h_minus_Mag_r_ij*h_minus_Mag_r_ij)/Mag_r_ij)*(r_ij);   // 1/mm^4 Vector

        /* Now apply the force to the two interacting bodies (Note the forces
        are equal and opposite). It should be noted that we don't actually
        apply the force to body B's jth particle. The reason for this is
        race conditions. Two threads can attempt to write
        the contact force to the same particle in body B at the same time.
        This causes a race condition (when run in parallel). This can be
        fixed using a critical region, but that's slow. To fix this, each
        thread gets its own 'Local contact force array'. each thread stores
        the force that it would apply to the particles in body B in this array.
        Once we have cycled through all of body A's particles, we then 'reduce'
        these force arrays. One by one, the threads add their contributions to
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
        F_Friction = ((-1*Simulation::Friction_Coefficient*Mag_F_Contact) / Relative_Velocity.Magnitude())*(Relative_Velocity);
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
  if(Contact_Flag == true) {
    for(i = 0; i < Num_Particles_B; i++) {
      Body_B[i].Force_Contact += Body_B_F_Contact_Local[i];                    //        : N Vector
      Body_B[i].Force_Friction += Body_B_F_Friction_Local[i];                  //        : N Vector
    } // for(i = 0; i < Num_Particles_B; i++) {
  } // if(Contact_Flag == true) {

  // Now free any dynamically allocated memory.
  delete [] Body_B_x;
  delete [] Body_B_F_Contact_Local;
  delete [] Body_B_F_Friction_Local;

  /* We put a barrier here so that each particle's contact and friction forces
  have been set before returning. The next function in the simulaiton (update_x)
  will only work correctly if each particle's contact and friction forces have
  been set. */
  #pragma omp barrier
} // void Body::Contact(Body & Body_A, Body & Body_B) {

struct Contact_Particle_Bucket {
  Array<unsigned> Array_A{};
  Array<unsigned> Array_B{};

  List<unsigned> List_A{};
  List<unsigned> List_B{};
}; // typdef struct Particle_Bucket {


void Body::Contact_New(Body & Body_A, Body & Body_B) {
  /* This function calculates the contact force between the particles in both A
  and those in body B.

  First, we partition the spatial domain. To do this, we first determine the
  maximum and minimum x, y, and z coordinates of the live (damage < 1) particles
  This gives us a cuboid in which the particles of the two bodies live. We
  then divide the x, y, and z dimensions of this cuboid into smaller cuboids,
  each one of which has a side length that is barely greater than the contact
  distance. We then allocate an array of buckets with one bucket per sub-cuboid.
  We then cycle through the particles of A and B, sorting them into their
  corresponding buckets (specifically the lists in those buckets).

  Once every particle has been sorted, the bucket lists are used to generate
  arrays in each bucket (we use lists at first because we don't know ahead of
  time how many particles will go into each bucket. We convert these lists into
  arrays so that we can access them quickly).

  We then apply the contact algorithm. The way that these buckets are
  set up, a particle can only contact particles that are in its bucket or a
  bucket which is adjacent to its bucket (27 total buckets to check). This
  greatly reduces the number of particles that we need to check against, thereby
  eliminating needless calculations. */

  //////////////////////////////////////////////////////////////////////////////
  /* Determine how many buckets we need */

  // Get paramaters from Body_A, Body_B
  const unsigned Num_Particles_A = Body_A.Get_Num_Particles();
  const unsigned Num_Particles_B = Body_B.Get_Num_Particles();

  // Set up variables.
  double x_max, x_min;
  double y_max, y_min;
  double z_max, z_min;

  Vector x = Body_A.Particles[0].Get_x();
  x_max = x[0];
  x_min = x_max;
  y_max = x[1];
  y_min = y_max;
  z_max = x[2];
  z_min = z_max;

  /* Determine the maximum and minimum x, y, z coordinate of the particles
  in the two bodies. */
  for(unsigned i = 0; i < Num_Particles_A; i++) {
    x = Body_A.Particles[i].Get_x();

    if(x[0] > x_max) { x_max = x[0]; }
    if(x[0] < x_min) { x_min = x[0]; }

    if(x[1] > y_max) { y_max = x[1]; }
    if(x[1] < y_min) { y_min = x[1]; }

    if(x[2] > z_max) { z_max = x[2]; }
    if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  for(unsigned i = 0; i < Num_Particles_B; i++) {
    x = Body_B.Particles[i].Get_x();

    if(x[0] > x_max) { x_max = x[0]; }
    if(x[0] < x_min) { x_min = x[0]; }

    if(x[1] > y_max) { y_max = x[1]; }
    if(x[1] < y_min) { y_min = x[1]; }

    if(x[2] > z_max) { z_max = x[2]; }
    if(x[2] < z_min) { z_min = x[2]; }
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {

  /* We want the dimension (in all three coordinate directions) of the
  sub-cuboids to be >= the contact distance. By doing this, particles
  can only compe into contact with particles in their bucket or in buckets
  that are adjacent to their bucket. In general, we want the buckets to be as
  small as possible (so that there are as few particles to check for contact as
  possible). Let's focus on the x coordinate. Let Nx denote the number of
  sub-cuboids in the x direction. We want Nx to be the largest natural number
  such that Contact_Distance <= (x_max - x_min)/Nx. A little though reveals
  that this occurs precisely when Nx = floor((x_max - x_min)/Contact_Distance).
  A similar result holds for the y and z directions. */
  const unsigned Nx = floor((x_max - x_min)/Simulation::Contact_Distance);
  const unsigned Ny = floor((y_max - y_min)/Simulation::Contact_Distance);
  const unsigned Nz = floor((z_max - z_min)/Simulation::Contact_Distance);

  /* Now that we know the number of buckets in each direction, we can allocate
  the buckets array! */
  Array<struct Contact_Particle_Bucket> Buckets{Nx*Ny*Nz};



  //////////////////////////////////////////////////////////////////////////////
  /* Now, sort the particles of A and B into their corresponding buckets. */

  /* First, we calculate some variables (see the next comment for an
  explanation) */
  const double sub_cuboid_x_dim = (x_max - x_min)/Nx;
  const double sub_cuboid_y_dim = (y_max - y_min)/Ny;
  const double sub_cuboid_z_dim = (y_max - z_min)/Nz;

  for(unsigned i = 0; i < Num_Particles_A; i++) {
    /* First, we need to determine which bucket our particle belongs in. Let's
    focus on the x coordinate. Each bucket has an x-dimension length of
    (x_max - x_min)/Nx, which we call sub_cuboid_x_dim. The x coordinates of
    the bucket for the ith particle is the number of units of length
    sub_cuboid_x_dim that fit between x_min and the x coordinate of the particle
    (think about it). A little thought reveals that this is precisely
    floor((Particle[i].x[0] - x_min)/((x_max - x_min)/Nx))
    A similar result holds for y and z. */
    x = Body_A.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/sub_cuboid_x_dim);
    unsigned ny = floor((x[1] - y_min)/sub_cuboid_y_dim);
    unsigned nz = floor((x[2] - z_min)/sub_cuboid_z_dim);

    Buckets[nx + ny*Nx + nz*Nx*Ny].List_A.Push_Front(i);
  } // for(unsigned i = 0; i < Num_Particles_A; i++) {

  for(unsigned i = 0; i < Num_Particles_B; i++) {
    x = Body_B.Particles[i].Get_x();
    unsigned nx = floor((x[0] - x_min)/sub_cuboid_x_dim);
    unsigned ny = floor((x[1] - y_min)/sub_cuboid_y_dim);
    unsigned nz = floor((x[2] - z_min)/sub_cuboid_z_dim);

    Buckets[nx + ny*Nx + nz*Nx*Ny].List_B.Push_Front(i);
  } // for(unsigned i = 0; i < Num_Particles_B; i++) {
} // void Body::Contact_New(Body & Body_A, Body & Body_B) {
