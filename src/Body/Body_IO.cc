#include "Body.h"
#include "Vector/Vector.h"
#include "Particle/Particle.h"
#include "Errors.h"
#include <stdio.h>
#include <cstring>
#include <string>

// Static function Prototypes
static void Add_Point_Data(FILE * File,
                           char * Weight_Name,
                           unsigned Num_Particles,
                           double * Data);

void Body::Print_Parameters(void) const {
  /* This function prints some of the key parameters about this body.

  This function is called by Simulation::Setup to check that everything ran
  correctly. */

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
    throw Not_A_Box(Buf);
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
  and the particle's mass, which is in units of grams. We need to divide the
  acceleration by 1,000,000 to make the units work. */
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
      Force_x_plus  += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
    } // for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {

  // -x face (i = 0)
  i = 0;
  for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_x_minus += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {

  // +y face (j = y_Side_len-1)
  j = (*this).Y_SIDE_LENGTH-1;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_y_plus  += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
    } //for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // -y face (j = 0)
  j = 0;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_y_minus += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
    } // for(k = 0; k < (*this).Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // +z face (k = Z_SIDE_LENGTH-1)
  k = (*this).Z_SIDE_LENGTH-1;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_z_plus  += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
    } // for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {

  // -z face (k = 0)
  k = 0;
  for(i = 0; i < (*this).X_SIDE_LENGTH; i++) {
    for(j = 0; j < (*this).Y_SIDE_LENGTH; j++) {
      unsigned index = i*(*this).Y_SIDE_LENGTH*(*this).Z_SIDE_LENGTH + k*(*this).Y_SIDE_LENGTH + j;
      Force_z_minus += Particles[index].Get_Mass()*(Particles[index].Get_a()/1000000.);
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
  /* This function prints (to a file) the position of every particle in (*this).
  This is used to visualize what the simulation is doing.

  This file SHOULD NOT be called by multiple threads. It is NOT thread safe. */

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
  /* Calculate particle data
  Here we calculate the particle data that the user has requested we print. The
  quantities that we can calculate point data for are:
    D:             Particle Damage
    Stretch_Max:   Maximum stretch (maximum eigenvalue of C)
    J:             det(F), local volume change ratio.
    F:             Deformation gradient.
    C:             Right Cauchy-Green Strain Tensor.
    E:             Green Strain tensor.
    P:             First Piola-Kirchoff Stress Tensor.
    T:             Cauchy Stress tensor.

  C, E, P, and T are all symmetric tensors. As such, we only print out their
  upper triangular components. F, by contrast, is not symmetric in general.
  Thus, we need to print out all 9 of its components.

  D, Stretch_Max, F, and P can acquired directly from the particle. The other
  quantities require calculations. C and E are functions of F. T is a function
  of T and P.

  For each quantity that the user wants, we save the relevant data into a set
  of dynamic arrays. Once this is finished, we write the components to the
  output file. */

  /* Create dynamic arrays for data that the user wants to print. */

  // Now print these values to the file.
  fprintf(File,"POINT_DATA %i\n", Num_Particles);
  char Weight_Name[15];

  // D
  if(Simulation::Print_Particle_D == true) {
    double * D = new double[Num_Particles];

    for(unsigned i = 0; i < Num_Particles; i++) { D[i] = Particles[i].Get_D(); }

    std::strcpy(Weight_Name, "D");
    Add_Point_Data(File, Weight_Name, Num_Particles, D);

    delete [] D;
  } // if(Simulation::Print_Particle_D == true) {


  // Stretch_Max
  if(Simulation::Print_Particle_Stretch_Max == true) {
    double * Stretch_Max = new double[Num_Particles];

    for(unsigned i = 0; i < Num_Particles; i++) { Stretch_Max[i] = Particles[i].Get_Stretch_M(); }

    std::strcpy(Weight_Name, "Stretch_Max");
    Add_Point_Data(File, Weight_Name, Num_Particles, Stretch_Max);

    delete [] Stretch_Max;
  } // if(Simulation::Print_Particle_Stretch_Max == true) {


  // J
  if(Simulation::Print_Particle_J == true) {
    double * J = new double[Num_Particles];

    Tensor F{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);
      J[i] = Determinant(F);
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "J");
    Add_Point_Data(File, Weight_Name, Num_Particles, J);
    delete [] J;
  } // if(Simulation::Print_Particle_J == true) {


  // F
  if(Simulation::Print_Particle_F == true) {
    double * F11 = new double[Num_Particles];
    double * F12 = new double[Num_Particles];
    double * F13 = new double[Num_Particles];

    double * F21 = new double[Num_Particles];
    double * F22 = new double[Num_Particles];
    double * F23 = new double[Num_Particles];

    double * F31 = new double[Num_Particles];
    double * F32 = new double[Num_Particles];
    double * F33 = new double[Num_Particles];

    Tensor F{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);

      F11[i] = F[0*3 + 0];
      F12[i] = F[0*3 + 1];
      F13[i] = F[0*3 + 2];

      F21[i] = F[1*3 + 0];
      F22[i] = F[1*3 + 1];
      F23[i] = F[1*3 + 2];

      F31[i] = F[2*3 + 0];
      F32[i] = F[2*3 + 1];
      F33[i] = F[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "F11");   Add_Point_Data(File, Weight_Name, Num_Particles, F11);
    std::strcpy(Weight_Name, "F12");   Add_Point_Data(File, Weight_Name, Num_Particles, F12);
    std::strcpy(Weight_Name, "F13");   Add_Point_Data(File, Weight_Name, Num_Particles, F13);

    std::strcpy(Weight_Name, "F21");   Add_Point_Data(File, Weight_Name, Num_Particles, F21);
    std::strcpy(Weight_Name, "F22");   Add_Point_Data(File, Weight_Name, Num_Particles, F22);
    std::strcpy(Weight_Name, "F23");   Add_Point_Data(File, Weight_Name, Num_Particles, F23);

    std::strcpy(Weight_Name, "F31");   Add_Point_Data(File, Weight_Name, Num_Particles, F31);
    std::strcpy(Weight_Name, "F32");   Add_Point_Data(File, Weight_Name, Num_Particles, F32);
    std::strcpy(Weight_Name, "F33");   Add_Point_Data(File, Weight_Name, Num_Particles, F33);

    delete [] F11;
    delete [] F12;
    delete [] F13;

    delete [] F21;
    delete [] F22;
    delete [] F23;

    delete [] F31;
    delete [] F32;
    delete [] F33;
  } // if(Simulation::Print_Particle_F == true) {


  // C
  if(Simulation::Print_Particle_C == true) {
    double * C11 = new double[Num_Particles];
    double * C12 = new double[Num_Particles];
    double * C13 = new double[Num_Particles];

    double * C22 = new double[Num_Particles];
    double * C23 = new double[Num_Particles];

    double * C33 = new double[Num_Particles];

    Tensor F{}, C{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);
      C = (F^T)*F;

      C11[i] = C[0*3 + 0];
      C12[i] = C[0*3 + 1];
      C13[i] = C[0*3 + 2];

      C22[i] = C[1*3 + 1];
      C23[i] = C[1*3 + 2];

      C33[i] = C[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "C11");   Add_Point_Data(File, Weight_Name, Num_Particles, C11);
    std::strcpy(Weight_Name, "C12");   Add_Point_Data(File, Weight_Name, Num_Particles, C12);
    std::strcpy(Weight_Name, "C13");   Add_Point_Data(File, Weight_Name, Num_Particles, C13);

    std::strcpy(Weight_Name, "C22");   Add_Point_Data(File, Weight_Name, Num_Particles, C22);
    std::strcpy(Weight_Name, "C23");   Add_Point_Data(File, Weight_Name, Num_Particles, C23);

    std::strcpy(Weight_Name, "C33");   Add_Point_Data(File, Weight_Name, Num_Particles, C33);

    delete [] C11;
    delete [] C12;
    delete [] C13;

    delete [] C22;
    delete [] C23;

    delete [] C33;
  } // if(Simulation::Print_Particle_C == true) {


  // B
  if(Simulation::Print_Particle_B == true) {
    double * B11 = new double[Num_Particles];
    double * B12 = new double[Num_Particles];
    double * B13 = new double[Num_Particles];

    double * B22 = new double[Num_Particles];
    double * B23 = new double[Num_Particles];

    double * B33 = new double[Num_Particles];

    Tensor F{}, B{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);
      B = F*(F^T);

      B11[i] = B[0*3 + 0];
      B12[i] = B[0*3 + 1];
      B13[i] = B[0*3 + 2];

      B22[i] = B[1*3 + 1];
      B23[i] = B[1*3 + 2];

      B33[i] = B[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "B11");   Add_Point_Data(File, Weight_Name, Num_Particles, B11);
    std::strcpy(Weight_Name, "B12");   Add_Point_Data(File, Weight_Name, Num_Particles, B12);
    std::strcpy(Weight_Name, "B13");   Add_Point_Data(File, Weight_Name, Num_Particles, B13);

    std::strcpy(Weight_Name, "B22");   Add_Point_Data(File, Weight_Name, Num_Particles, B22);
    std::strcpy(Weight_Name, "B23");   Add_Point_Data(File, Weight_Name, Num_Particles, B23);

    std::strcpy(Weight_Name, "B33");   Add_Point_Data(File, Weight_Name, Num_Particles, B33);

    delete [] B11;
    delete [] B12;
    delete [] B13;

    delete [] B22;
    delete [] B23;

    delete [] B33;
  } // if(Simulation::Print_Particle_B == true) {


  // E
  if(Simulation::Print_Particle_E == true) {
    double * E11 = new double[Num_Particles];
    double * E12 = new double[Num_Particles];
    double * E13 = new double[Num_Particles];

    double * E22 = new double[Num_Particles];
    double * E23 = new double[Num_Particles];

    double * E33 = new double[Num_Particles];

    Tensor F{}, E{};
    Tensor I{1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);
      E = (1./2.)*((F^T)*F - I);

      E11[i] = E[0*3 + 0];
      E12[i] = E[0*3 + 1];
      E13[i] = E[0*3 + 2];

      E22[i] = E[1*3 + 1];
      E23[i] = E[1*3 + 2];

      E33[i] = E[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "E11");   Add_Point_Data(File, Weight_Name, Num_Particles, E11);
    std::strcpy(Weight_Name, "E12");   Add_Point_Data(File, Weight_Name, Num_Particles, E12);
    std::strcpy(Weight_Name, "E13");   Add_Point_Data(File, Weight_Name, Num_Particles, E13);

    std::strcpy(Weight_Name, "E22");   Add_Point_Data(File, Weight_Name, Num_Particles, E22);
    std::strcpy(Weight_Name, "E23");   Add_Point_Data(File, Weight_Name, Num_Particles, E23);

    std::strcpy(Weight_Name, "E33");   Add_Point_Data(File, Weight_Name, Num_Particles, E33);

    delete [] E11;
    delete [] E12;
    delete [] E13;

    delete [] E22;
    delete [] E23;

    delete [] E33;
  } // if(Simulation::Print_Particle_E == true) {


  // P
  if(Simulation::Print_Particle_P == true) {
    double * P11 = new double[Num_Particles];
    double * P12 = new double[Num_Particles];
    double * P13 = new double[Num_Particles];

    double * P22 = new double[Num_Particles];
    double * P23 = new double[Num_Particles];

    double * P33 = new double[Num_Particles];

    Tensor P{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      P = Particles[i].Get_P();

      P11[i] = P[0*3 + 0];
      P12[i] = P[0*3 + 1];
      P13[i] = P[0*3 + 2];

      P22[i] = P[1*3 + 1];
      P23[i] = P[1*3 + 2];

      P33[i] = P[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "P11");   Add_Point_Data(File, Weight_Name, Num_Particles, P11);
    std::strcpy(Weight_Name, "P12");   Add_Point_Data(File, Weight_Name, Num_Particles, P12);
    std::strcpy(Weight_Name, "P13");   Add_Point_Data(File, Weight_Name, Num_Particles, P13);

    std::strcpy(Weight_Name, "P22");   Add_Point_Data(File, Weight_Name, Num_Particles, P22);
    std::strcpy(Weight_Name, "P23");   Add_Point_Data(File, Weight_Name, Num_Particles, P23);

    std::strcpy(Weight_Name, "P33");   Add_Point_Data(File, Weight_Name, Num_Particles, P33);

    delete [] P11;
    delete [] P12;
    delete [] P13;

    delete [] P22;
    delete [] P23;

    delete [] P33;
  } // if(Simulation::Print_Particle_P == true) {


  // T
  if(Simulation::Print_Particle_T == true) {
    double * T11 = new double[Num_Particles];
    double * T12 = new double[Num_Particles];
    double * T13 = new double[Num_Particles];

    double * T22 = new double[Num_Particles];
    double * T23 = new double[Num_Particles];

    double * T33 = new double[Num_Particles];

    double J;
    Tensor F{}, P{}, T{};
    for(unsigned i = 0; i < Num_Particles; i++) {
      F = Particles[i].Get_F((*this).F_Index);
      P = Particles[i].Get_P();
      J = Determinant(F);
      T = (P*Transpose(F))/J;

      T11[i] = T[0*3 + 0];
      T12[i] = T[0*3 + 1];
      T13[i] = T[0*3 + 2];

      T22[i] = T[1*3 + 1];
      T23[i] = T[1*3 + 2];

      T33[i] = T[2*3 + 2];
    } // for(unsigned i = 0; i < Num_Particles; i++) {

    std::strcpy(Weight_Name, "T11");   Add_Point_Data(File, Weight_Name, Num_Particles, T11);
    std::strcpy(Weight_Name, "T12");   Add_Point_Data(File, Weight_Name, Num_Particles, T12);
    std::strcpy(Weight_Name, "T13");   Add_Point_Data(File, Weight_Name, Num_Particles, T13);

    std::strcpy(Weight_Name, "T22");   Add_Point_Data(File, Weight_Name, Num_Particles, T22);
    std::strcpy(Weight_Name, "T23");   Add_Point_Data(File, Weight_Name, Num_Particles, T23);

    std::strcpy(Weight_Name, "T33");   Add_Point_Data(File, Weight_Name, Num_Particles, T33);

    delete [] T11;
    delete [] T12;
    delete [] T13;

    delete [] T22;
    delete [] T23;

    delete [] T33;
  } // if(Simulation::Print_Particle_T == true) {


  // Free the file
  fclose(File);
} // void Body::Export_Particle_Positions(void) {



static void Add_Point_Data(FILE * File, char * Weight_Name, unsigned Num_Particles, double * Data) {
  /* This function adds some data to a vtk type file. This function exists
  because VTK expects a particular format, and I thought it'd be better to write
  one function that does that rather than to have a bunch. */

  // Print header.
  fprintf(File, "SCALARS ");
  fprintf(File, Weight_Name);
  fprintf(File, " float\n");
  fprintf(File, "LOOKUP_TABLE default\n");

  // Now print supplied data to file
  for(unsigned i = 0; i < Num_Particles; i++) {
    fprintf(File,"\t %8.3f\n", Data[i]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {
} // static void Add_Point_Data(FILE * File, char * Weight_Name, unsigned Num_Particles, double * Data) const {
