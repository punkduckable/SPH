#include "Body.h"
#include "Vector/Vector.h"
#include "Particle/Particle.h"
#include "Errors.h"
#include <stdio.h>
#include <cstring>
#include <string>

void Body::Print_Parameters(void) const {
  printf(         "Name:                         %s\n",    Name.c_str());
  printf(         "Is a Box:                     %u\n",    (unsigned)Is_Box);
  if(Is_Box == true) {
    printf(       "X side length:                %u\n",    X_SIDE_LENGTH);
    printf(       "Y side length:                %u\n",    Y_SIDE_LENGTH);
    printf(       "Z side length:                %u\n",    Z_SIDE_LENGTH);
  } // if(Is_Box) {
  printf(         "Number of particles:          %u\n",    Num_Particles);
  printf(         "Particles address:            %p\n",    Particles);
  printf(         "Inter particle spacing (mm):  %lf\n",   Inter_Particle_Spacing);
  printf(         "Support Radius (mm):          %lf\n",   Support_Radius);
  printf(         "Shape Function Amplitude:     %lf\n",   Shape_Function_Amplitude);

  // Material parameters
  printf(         "Lame (Mpa):                   %lf\n",   Body_Material.Lame);
  printf(         "mu0 (Shear modulus) (Mpa:     %lf\n",   Body_Material.mu0);
  printf(         "E (Young's modulus) (Mpa):    %lf\n",   Body_Material.E);
  printf(         "Density (g/mm^3):             %lf\n",   Body_Material.E);

  // Other
  printf(         "Gravity Enabled:              %u\n",    Gravity_Enabled);
  printf(         "mu (Viscosity):               %lf\n",   mu);

  printf(         "Is damageable:                %u\n",    Is_Damageable);
  if(Is_Damageable == true) {
    printf(       "Tau (Damage rate):            %lf\n\n", Tau);
  } // if(Is_Damageable == true) {

} // void Body::Print_Parameters(void) const {



void Body::Export_Body_Forces(const unsigned time_steps) {
  /* This function is used to find and print the forces applied to a body.

  This function can NOT be called by multiple threads at once (this
  function is not thread safe). */

  #if defined(IO_MONITOR)
    printf("Exporting Forces for %s\n",(*this).Name.c_str());
  #endif

  // First, open the file.
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name.c_str();
  File_Path +=  "_Forces.txt";

  FILE * File;
  if(Times_Printed_Body_Forces == 0) { File = fopen(File_Path.c_str(),"w"); }
  else {                               File = fopen(File_Path.c_str(),"a"); }

  // Make sure we could open the file.
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Cant Open File Exception: Thrown by Body::Export_Body_Forces\n"
            "For some reason, ./IO/Force_Files/%s_Forces.txt wouldn't open :(\n",
            (*this).Name.c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {

  /* Calculate the Elastic, Viscosity, Contact, Friction, and Hourglass forces
  acting acting on the body. To do this, we add up the corresponding forces in
  each particle in the body.

  Note: the internal force combines the force due to the elastic and viscous
  parts of the constitutive law. I want to look at the contribuitions from
  both pieces. Thankfully, we calculate the viscous contributions on their own.
  Thus, to get the elastic terms, I just subtract the viscosity force from the
  internal one.

  Note: we count damaged particles in this calculation. Thus, this
  calculation may give inaccurate results if some of the body's particles are
  damaged */
  Vector Internal_Force  = {0, 0, 0};
  Vector Elastic_Force   = {0, 0, 0};
  Vector Viscosity_Force = {0, 0, 0};
  Vector Contact_Force   = {0, 0, 0};
  Vector Friction_Force  = {0, 0, 0};
  Vector Hourglass_Force = {0, 0, 0};
  Vector Net_Force       = {0, 0, 0};

  for(unsigned i = 0; i < Num_Particles; i++) {
    /*  Note: to calculate the total force, we use a, which is in units of mm/s^2,
    and the particle's mass, which is in units of grams. We need some conversion
    factors to make this work. */
    Internal_Force  += Particles[i].Get_Force_Internal();
    Viscosity_Force += Particles[i].Get_Force_Viscosity();
    Contact_Force   += Particles[i].Get_Force_Contact();
    Friction_Force  += Particles[i].Get_Force_Friction();
    Hourglass_Force += Particles[i].Get_Force_Hourglass();
    Net_Force       += (Particles[i].Get_Mass()/1000.)*(Particles[i].Get_a()/1000.);
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  Elastic_Force = Internal_Force - Viscosity_Force;

  /* Print the results to file. If we're on the first time step, then we need
  to print a header. Otherwise, just print the forces! */
  if(Times_Printed_Body_Forces == 0) {
    fprintf(File,"Time Steps |");
    fprintf(File,"            Elastic Force (N)             |");
    fprintf(File,"            Viscous Force (N)             |");
    fprintf(File,"            Contact Force (N)             |");
    fprintf(File,"            Friction Force (N)            |");
    fprintf(File,"            Hourglass Force (N)           |");
    fprintf(File,"            Net Force (N)                 \n");
  } // if(Times_Printed_Body_Forces == 0) {

  fprintf(File,"%10d | ", time_steps);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Elastic_Force[0],   Elastic_Force[1],   Elastic_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Viscosity_Force[0], Viscosity_Force[1], Viscosity_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Contact_Force[0],   Contact_Force[1],   Contact_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Friction_Force[0],  Friction_Force[1],  Friction_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Hourglass_Force[0], Hourglass_Force[1], Hourglass_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e>\n",   Net_Force[0],       Net_Force[1],       Net_Force[2]);

  // Now close the file.
  fclose(File);

  // Increment the number of times that we're printed Body force data.
  (*this).Times_Printed_Body_Forces++;
} // void Body::Export_Body_Forces(const unsigned time_steps) {



void Body::Export_Body_Torques(const unsigned time_steps) {
  /* This function is used to find and print the torques applied to a body.

  This function can NOT be called by multiple threads at once (this
  function is not thread safe). */

  #if defined(IO_MONITOR)
    printf("Exporting Torques for %s\n",(*this).Name.c_str());
  #endif

  // First, open the file.
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name.c_str();
  File_Path +=  "_Torques.txt";

  FILE * File;
  if(Times_Printed_Body_Torques == 0) { File = fopen(File_Path.c_str(),"w"); }
  else {                                File = fopen(File_Path.c_str(),"a"); }

  // Make sure we could open the file.
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Cant Open File Exception: Thrown by Body::Export_Body_Torques\n"
            "For some reason, ./IO/Force_Files/%s_Torques.txt wouldn't open :(\n",
            (*this).Name.c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {

  /* Next, we need to find the centroid of the body. Since every particle in
  the body has the same mass, this is equivalent to just finding the average
  position of the particles in the body.

  Note: we still count damaged particles in this calculation. Thus, this
  calculation may give inaccurate results if some of the body's particles are
  damaged */
  Vector Centroid = {0, 0, 0};
  for(unsigned i = 0; i < Num_Particles; i++) {
    Centroid += (*this).Particles[i].Get_x();
  } // for(unsigned i = 0; i < Num_Particles; i++) {
  Centroid = Centroid/Num_Particles;

  /* Calculate the Torques due to Elastic, Viscosity, Contact, Friction, and
  Hourglass forces. To do this, we add up the corresponding Torques in
  each particle in the body. */
  Vector Internal_Torque  = {0, 0, 0};
  Vector Elastic_Torque   = {0, 0, 0};
  Vector Viscosity_Torque = {0, 0, 0};
  Vector Contact_Torque   = {0, 0, 0};
  Vector Friction_Torque  = {0, 0, 0};
  Vector Hourglass_Torque = {0, 0, 0};
  Vector Net_Torque       = {0, 0, 0};

  for(unsigned i = 0; i < Num_Particles; i++) {
    /* First, calculate the displacement between the Body's centroid and the
    ith particle's spatial position. */
    Vector Displacement = (*this).Particles[i].Get_x() - Centroid;

    /* Next calculate the torque due to the ith particle from each force acting
    on the particle. In general, the torque on the body due to a force F acting
    on a point with some displacement from the body's centroid is
    F x Displacement.

    Note: to calculate the total force, we use a, which is in units of mm/s^2,
    and the particle's mass, which is in units of grams. We need some conversion
    factors to make this work. */
    Internal_Torque  += Cross_Product(Displacement, Particles[i].Get_Force_Internal());
    Viscosity_Torque += Cross_Product(Displacement, Particles[i].Get_Force_Viscosity());
    Contact_Torque   += Cross_Product(Displacement, Particles[i].Get_Force_Contact());
    Friction_Torque  += Cross_Product(Displacement, Particles[i].Get_Force_Friction());
    Hourglass_Torque += Cross_Product(Displacement, Particles[i].Get_Force_Hourglass());
    Net_Torque       += Cross_Product(Displacement, (Particles[i].Get_Mass()/1000.)*(Particles[i].Get_a()/1000.));
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  Elastic_Torque = Internal_Torque - Viscosity_Torque;

  /* Print the results to file. If we're on the first time step, then we need
  to print a header. Otherwise, just print the torques! */
  if(Times_Printed_Body_Torques == 0) {
    fprintf(File,"Time Steps |");
    fprintf(File,"           Elastic Torque (N)             |");
    fprintf(File,"           Viscous Torque (N)             |");
    fprintf(File,"           Contact Torque (N)             |");
    fprintf(File,"           Friction Torque (N)            |");
    fprintf(File,"           Hourglass Torque (N)           |");
    fprintf(File,"           Net Torque (N)                 \n");
  } //   if(Times_Printed_Body_Torques == 0) {

  fprintf(File,"%10d | ", time_steps);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Elastic_Torque[0],   Elastic_Torque[1],   Elastic_Torque[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Viscosity_Torque[0], Viscosity_Torque[1], Viscosity_Torque[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Contact_Torque[0],   Contact_Torque[1],   Contact_Torque[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Friction_Torque[0],  Friction_Torque[1],  Friction_Torque[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Hourglass_Torque[0], Hourglass_Torque[1], Hourglass_Torque[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e>\n",   Net_Torque[0],       Net_Torque[1],       Net_Torque[2]);

  // Now close the file.
  fclose(File);

  // Increment the number of times that we're printed Body force data.
  (*this).Times_Printed_Body_Torques++;
} // void Body::Export_Body_Torques(const unsigned time_steps) {



void Body::Export_Box_Boundary_Forces(const unsigned time_steps) {
  /* This function is used to print the resultant forces applied to each of
  a box's six boundaries. Since this information only really makes sense for
  boxes, this function is only available to boxes. */

  #if defined(IO_MONITOR)
    printf("Exporting Boundary Forces for %s\n",(*this).Name.c_str());
  #endif

  // Verify that this is a box!
  if((*this).Is_Box == false) {
    char Buf[500];
    sprintf(Buf,
            "Not A Box Exception: thrown by Body::Export_Box_Boundary_Forces\n"
            "Body %s tried to use this function, but %s is not a box! This function\n"
            "can only be called by boxes!\n",
            (*this).Name.c_str(), (*this).Name.c_str());
    throw Cant_Open_File(Buf);
  } // if((*this).Is_Box == false) {

  // First, open the file.
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name.c_str();
  File_Path +=  "_Boundary_Forces.txt";

  FILE * File;
  if(Times_Printed_Box_Boundary_Forces == 0) { File = fopen(File_Path.c_str(),"w"); }
  else {                                       File = fopen(File_Path.c_str(),"a"); }

  // Make sure we could open the file.
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Cant Open File Exception: Thrown by Body::Export_Box_Boundary_Forces\n"
            "For some reason, ./IO/Force_Files/%s_Boundary_Forces.txt wouldn't open :(\n",
            (*this).Name.c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {

  /* Next, we need to find the net forces applied to each of the six boundaries
  of the box. To do this, we add up the net force applied to each particle in
  each of the 6 boundaries.

  Note: to calculate the total force, we use a, which is in units of mm/s^2,
  and the particle's mass, which is in units of grams. We need some conversion
  factors to make this work. */
  Vector Force_x_plus  = {0, 0, 0};
  Vector Force_x_minus = {0, 0, 0};
  Vector Force_y_plus  = {0, 0, 0};
  Vector Force_y_minus = {0, 0, 0};
  Vector Force_z_plus  = {0, 0, 0};
  Vector Force_z_minus = {0, 0, 0};

  unsigned i,j,k;

  // +x face (i = X_SIDE_LENGTH-1)
  i = (*this).X_SIDE_LENGTH-1;
  for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_x_plus  += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);
    } // for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {

  // -x face (i = 0)
  i = 0;
  for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_x_minus += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {

  // +y face (j = y_Side_len-1)
  j = (*this).Y_SIDE_LENGTH-1;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_y_plus  += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);
    } //for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // -y face (j = 0)
  j = 0;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_y_minus += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);
    } // for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // +z face (k = Z_SIDE_LENGTH-1)
  k = (*this).Z_SIDE_LENGTH-1;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_z_plus  += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);
    } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // -z face (k = 0)
  k = 0;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_z_minus += (Particles[index].Get_Mass()/1000.)*(Particles[index].Get_a()/1000.);
    } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  /* Print the results to file. If we're on the first time step, then we need
  to print a header. Otherwise, just print the forces! */
  if(Times_Printed_Box_Boundary_Forces == 0) {
    fprintf(File,"Time Steps |");
    fprintf(File,"          x+ Boundary Force (N)           |");
    fprintf(File,"          x- Boundary Force (N)           |");
    fprintf(File,"          y+ Boundary Force (N)           |");
    fprintf(File,"          y- Boundary Force (N)           |");
    fprintf(File,"          z+ Boundary Force (N)           |");
    fprintf(File,"          z- Boundary Force (N)          \n");
  } //   if(Times_Printed_Box_Boundary_Forces == 0) {

  fprintf(File,"%10d | ", time_steps);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Force_x_plus[0],  Force_x_plus[1],  Force_x_plus[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Force_x_minus[0], Force_x_minus[1], Force_x_minus[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Force_y_plus[0],  Force_y_plus[1],  Force_y_plus[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Force_y_minus[0], Force_y_minus[1], Force_y_minus[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Force_z_plus[0],  Force_z_plus[1],  Force_z_plus[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e>\n",   Force_z_minus[0], Force_z_minus[1], Force_z_minus[2]);

  // Now close the file.
  fclose(File);

  // Increment the number of times that we're printed Body force data.
  (*this).Times_Printed_Box_Boundary_Forces++;
} // void Body::Export_Box_Boundary_Forces(const unsigned time_steps) {



void Body::Export_Particle_Forces(void) {
  /* This function is used to print the forces applied to each particle in a
  body.

  This function can NOT be called by multiple threads at once (this
  function is not thread safe). */

  #if defined(IO_MONITOR)
    printf("Exporting particle forces for %s\n",(*this).Name.c_str());
  #endif

  // Create a file path for the new file (based on the Body's name
  // and time_step)
  char Buf[10];
  sprintf(Buf,"%05u.txt",Times_Printed_Particle_Forces);
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name;
  File_Path +=  "_Force_";
  File_Path +=  Buf;

  // Now open the file.
  FILE * File = fopen(File_Path.c_str(), "w");
  if(File == nullptr) {
    char Error_Buf[500];
    sprintf(Error_Buf,
            "Cant Open File Exception: Thrown by Body::Export_Particle_Forces\n"
            "For some reason, ./IO/Force_Files/%s_Force_%s won't open :(\n",
            (*this).Name.c_str(),
            Buf);
    throw Cant_Open_File(Error_Buf);
  } // if(File == nullptr) {

  // Increment the number of times that we're printed particle force data.
  (*this).Times_Printed_Particle_Forces++;

  // Print header.
  fprintf(File,"  ID  |");
  fprintf(File," Particle Pos  |");
  fprintf(File,"        Internal Force        |");
  fprintf(File,"        Viscous Force         |");
  fprintf(File,"        Contact Force         |");
  fprintf(File,"        Friction Force        |");
  fprintf(File,"        Hourglass Force       |");
  fprintf(File,"\n");

  // Cycle through particles, print spacial positions, forces for each particle
  for(unsigned i = 0; i < Num_Particles; i++) {
    fprintf(File,"%6u|", Particles[i].Get_ID());
    fprintf(File,"%4.1f,%4.1f,%4.1f | ",    Particles[i].X[0],               Particles[i].X[1],               Particles[i].X[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Internal[0],  Particles[i].Force_Internal[1],  Particles[i].Force_Internal[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Viscosity[0], Particles[i].Force_Viscosity[1], Particles[i].Force_Viscosity[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Contact[0],   Particles[i].Force_Contact[1],   Particles[i].Force_Contact[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Friction[0],  Particles[i].Force_Friction[1],  Particles[i].Force_Friction[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e>\n",   Particles[i].Force_Hourglass[0], Particles[i].Force_Hourglass[1], Particles[i].Force_Hourglass[2]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Body::Export_Particle_Forces(void) {



void Body::Export_Particle_Positions(void) {
  #if defined(IO_MONITOR)
    printf("Exporting particle positions for %s\n",(*this).Name.c_str());
  #endif

  // Set up file
  char Buf[10];
  sprintf(Buf,"%05u.vtk",Times_Printed_Particle_Positions);
  std::string File_Path = "./IO/Position_Files/";
  File_Path += (*this).Name;
  File_Path += "_positions_";
  File_Path += Buf;
  FILE * File = fopen(File_Path.c_str(), "w");
  if(File == nullptr) {
    char Error_Buf[500];
    sprintf(Error_Buf,
            "Cant Open File Exception: Thrown by Body::Export_Particle_Positions\n"
            "For some reason, ./IO/Position_Files/%s_positions_%s won't open :(\n",
            (*this).Name.c_str(),
            Buf);
    throw Cant_Open_File(Error_Buf);
  } // if(File == nullptr) {

  // Increment the number of times that we're printed particle positio data.
  Times_Printed_Particle_Positions++;



  //////////////////////////////////////////////////////////////////////////////
  // Print file header
  fprintf(File,"%s\n","# vtk DataFile Version 3.0");
  fprintf(File,"%s\n","test_file");
  fprintf(File,"%s\n","ASCII");
  fprintf(File,"%s\n","DATASET POLYDATA");
  fprintf(File,"POINTS %i float\n",Num_Particles);



  //////////////////////////////////////////////////////////////////////////////
  // Cycle through particles, print spacial positions of each particle
  Vector x;
  for(unsigned i = 0; i < Num_Particles; i++) {
    x = Particles[i].Get_x();
    fprintf(File,"%8.3f \t %8.3f \t %8.3f\n",x[0], x[1], x[2]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {



  //////////////////////////////////////////////////////////////////////////////
  /* Find the components of S and E for each particle
  We Calculate each of these by first finding P and F for each particle. We then
  use these quantities to calculate Sigma and E. The components of each of these
  tensors are then stored in dynamic arrays. Once this is finished, we write the
  components to the output file.

  We choose to store the components in dynamic arrays before writing to the file
  to improve performnace. Each particle's P, F tensors are in distinct memory
  locations. If we were to collect each component at a time, we'd have to read
  in a new cache line for each component from each particle. In other words,
  we'd only get one double from each cache line. This is bad. To improve this,
  we only pull from the particles once. We pull in the entire tensor, then write
  its components to the dynamic arrays. We use all 6 of each tensor's components
  in each iteration. This means that we can make better usage of cache lines.

  Once the components have been written to the arrays, tensor components
  belonging to adjacent particles are next to each ther in memory. For example,
  in the P11 array, P11 component of the 100th particle is stored right next to
  the P11 component of the 101th particle. when we write to the file, we pull
  from the dynamic arrays. This allows us to use cache lines efficiently.

  This improves performance.
  */

  /* Create dynamic arrays for components of S, E (note, both are symmetric, so
  we only need to store 6 components) and J (det F)*/

  double * LamM = new double[Num_Particles];
  //double * LamH = new double[Num_Particles];
  //double * LamC = new double[Num_Particles];
  double * D = new double[Num_Particles];
  /*
  double * S11 = new double[Num_Particles];
  double * S22 = new double[Num_Particles];
  double * S33 = new double[Num_Particles];
  double * S21 = new double[Num_Particles];
  double * S31 = new double[Num_Particles];
  double * S32 = new double[Num_Particles];
  */
  //double * E11 = new double[Num_Particles];
  //double * E22 = new double[Num_Particles];
  //double * E33 = new double[Num_Particles];
  //double * E21 = new double[Num_Particles];
  //double * E31 = new double[Num_Particles];
  //double * E32 = new double[Num_Particles];

  //double * J = new double[Num_Particles];

  Tensor F, P, S, E;
  Tensor I{1,0,0,
           0,1,0,
           0,0,1};

  for(unsigned i = 0; i < Num_Particles; i++) {
    LamM[i] = Particles[i].Get_Stretch_M();
    //LamH[i] = Particles[i].Get_Stretch_H();
    //LamC[i] = Particles[i].Get_Stretch_Critical();
    D[i] = Particles[i].Get_D();

    // Get F, P from current particle
    //F = Particles[i].Get_F();
    //P = Particles[i].Get_P();

    // Use F to calculate determinant (J)
    //J[i] = Determinant(F);

    // Calculate S from P.
    //S = P*(F^(T))/J[i];

    // Get components of S
    /*
    S11[i] = S[3*0 + 0];
    S22[i] = S[3*1 + 1];
    S33[i] = S[3*2 + 2];
    S21[i] = S[3*1 + 0];
    S31[i] = S[3*2 + 0];
    S32[i] = S[3*2 + 1];
    */

    // Now calculate E = (1/2)(C-I)
    //E = (1./2.)*((F^T)*F - I);

    // Now get components of E
    //E11[i] = E[3*0 + 0];
    //E22[i] = E[3*1 + 1];
    //E33[i] = E[3*2 + 2];
    //E21[i] = E[3*1 + 0];
    //E31[i] = E[3*2 + 0];
    //E32[i] = E[3*2 + 1];
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  // Now print these values to the file.
  fprintf(File,"POINT_DATA %i\n", Num_Particles);
  char Weight_Name[5];

  /* Damage parameters  */

  std::strcpy(Weight_Name, "LamM");
  Add_Point_Data(File, Weight_Name, Num_Particles, LamM);

  //std::strcpy(Weight_Name, "LamH");
  //Add_Point_Data(File, Weight_Name, Num_Particles, LamH);

  //std::strcpy(Weight_Name, "LamC");
  //Add_Point_Data(File, Weight_Name, Num_Particles, LamC);

  std::strcpy(Weight_Name, "D");
  Add_Point_Data(File, Weight_Name, Num_Particles, D);
  /* Components of S */
  /*
  std::strcpy(Weight_Name, "S11");
  Add_Point_Data(File, Weight_Name, Num_Particles, S11);

  std::strcpy(Weight_Name, "S22");
  Add_Point_Data(File, Weight_Name, Num_Particles, S22);

  std::strcpy(Weight_Name, "S33");
  Add_Point_Data(File, Weight_Name, Num_Particles, S33);

  std::strcpy(Weight_Name, "S21");
  Add_Point_Data(File, Weight_Name, Num_Particles, S21);

  std::strcpy(Weight_Name, "S31");
  Add_Point_Data(File, Weight_Name, Num_Particles, S31);

  std::strcpy(Weight_Name, "S32");
  Add_Point_Data(File, Weight_Name, Num_Particles, S32);
  */

  /* Components of E*/
  //std::strcpy(Weight_Name, "E11");
  //Add_Point_Data(File, Weight_Name, Num_Particles, E11);

  //std::strcpy(Weight_Name, "E22");
  //Add_Point_Data(File, Weight_Name, Num_Particles, E22);
  /*
  std::strcpy(Weight_Name, "E33");
  Add_Point_Data(File, Weight_Name, Num_Particles, E33);

  std::strcpy(Weight_Name, "E21");
  Add_Point_Data(File, Weight_Name, Num_Particles, E21);

  std::strcpy(Weight_Name, "E31");
  Add_Point_Data(File, Weight_Name, Num_Particles, E31);

  std::strcpy(Weight_Name, "E32");
  Add_Point_Data(File, Weight_Name, Num_Particles, E32);
  */

  /* J */
  /*
  std::strcpy(Weight_Name, "J");
  Add_Point_Data(File, Weight_Name, Num_Particles, J);
  */

  // Deallocate dynamic arrays
  delete [] LamM;
  //delete [] LamH;
  //delete [] LamC;
  delete [] D;
  /*
  delete [] S11;
  delete [] S22;
  delete [] S33;
  delete [] S21;
  delete [] S31;
  delete [] S32;
  */
  //delete [] E11;
  //delete [] E22;
  //delete [] E33;
  //delete [] E21;
  //delete [] E31;
  //delete [] E32;

  //delete [] J;

  // Free the file
  fclose(File);
} // void Body::Export_Particle_Positions(void) {



void Body::Add_Point_Data(FILE * File, char * Weight_Name, unsigned Num_Particles, double * Data) const {
  // Print header.
  fprintf(File, "SCALARS ");
  fprintf(File, Weight_Name);
  fprintf(File, " float\n");
  fprintf(File, "LOOKUP_TABLE default\n");

  // Now print supplied data to file
  for(unsigned i = 0; i < Num_Particles; i++) {
    fprintf(File,"\t %8.3f\n", Data[i]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {
} // void Body::Add_Point_Data(FILE * File, char * Weight_Name, unsigned Num_Particles, double * Data) const {
