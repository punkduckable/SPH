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
  printf(         "Partciles Array address:      %p\n",    Particles);
  printf(         "Inter particle spacing:       %lf\n",   Inter_Particle_Spacing);
  printf(         "h:                            %lf\n",   h);
  printf(         "Support Radius:               %u\n",    Support_Radius);
  printf(         "Shape Function Amplitude:     %lf\n",   Shape_Function_Amplitude);
  printf(         "Lame:                         %lf\n",   Body_Material.Lame);
  printf(         "mu0 (Shear modulus):          %lf\n",   Body_Material.mu0);
  printf(         "Gravity Enabled:              %u\n",    Gravity_Enabled);
  printf(         "mu (Viscosity):               %lf\n",   mu);
  printf(         "E (Young's modulus):          %lf\n",   Body_Material.E);
  printf(         "Tau (Damage rate):            %lf\n\n", Tau);
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

  /* Calculate the Internal, Viscosity, Contact, Friction, and Hourglass forces
  acting acting on the body. To do this, we add up the corresponding forces in
  each particle in the body. */
  Vector Internal_Force  = {0, 0, 0};
  Vector Viscosity_Force = {0, 0, 0};
  Vector Contact_Force   = {0, 0, 0};
  Vector Friction_Force  = {0, 0, 0};
  Vector Hourglass_Force = {0, 0, 0};
  Vector Net_Force       = {0, 0, 0};

  for(unsigned i = 0; i < Num_Particles; i++) {
    Internal_Force  += Particles[i].Get_Force_Internal();
    Viscosity_Force += Particles[i].Get_Force_Viscosity();
    Contact_Force   += Particles[i].Get_Force_Contact();
    Friction_Force  += Particles[i].Get_Force_Friction();
    Hourglass_Force += Particles[i].Get_Force_Hourglass();
    Net_Force       += (Particles[i].Get_Mass()/1000.)*(Particles[i].Get_a()/1000.);     // Particle mass is in g, we want it in Kg. a is in mm/s^2, we want it in m/s^2.
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  /* Print the results to file. If we're on the first time step, then we need
  to print a header. Otherwise, just print the forces! */
  if(Times_Printed_Body_Forces == 0) {
    fprintf(File,"Time Steps |");
    fprintf(File,"            Internal Force (N)            |");
    fprintf(File,"            Viscous Force (N)             |");
    fprintf(File,"            Contact Force (N)             |");
    fprintf(File,"            Friction Force (N)            |");
    fprintf(File,"            Hourglass Force (N)           |");
    fprintf(File,"            Net Force (N)                 \n");
  } //   if(Times_Printed_Body_Forces == 0) {

  fprintf(File,"%10d | ", time_steps);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Internal_Force[0],  Internal_Force[1],  Internal_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Viscosity_Force[0], Viscosity_Force[1], Viscosity_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Contact_Force[0],   Contact_Force[1],   Contact_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Friction_Force[0],  Friction_Force[1],  Friction_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e> | ",  Hourglass_Force[0], Hourglass_Force[1], Hourglass_Force[2]);
  fprintf(File,"<%12.5e,%12.5e,%12.5e>\n",   Net_Force[0],       Net_Force[1],       Net_Force[2]);

  // Now close the file.
  fclose(File);

  // Increment the number of times that we're printed Body force data.
  Times_Printed_Body_Forces++;
} // void Body::Export_Body_Forces(const unsigned time_steps) {



void Body::Export_Particle_Forces(void) {
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
  Times_Printed_Particle_Forces++;

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
  } // for(unsinged int i = 0; i < Num_Particles; i++) {

  // Now print these values to the file.
  fprintf(File,"POINT_DATA %i\n", Num_Particles);
  char Weight_Name[5];

  /* Damage paramaters */

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
