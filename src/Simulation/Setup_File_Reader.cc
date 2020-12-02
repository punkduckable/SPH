#include "Simulation.h"
#include "IO/IO_Ops.h"
#include "Errors.h"
#include "Body/Body.h"
#include "Vector/Vector.h"
#include "Array.h"
#include <fstream>
#include <string>
#include <stdio.h>

// Prototypes for functions that are local to this file.
namespace Simulation {
  void Load_Body_From_Setup_File(Body & Body_In, const unsigned i);
  Vector Parse_BC(std::string & Line);
  void Print_BC(Vector & BC_In);
} // namespace Simulation {



Body* Simulation::Load_Setup_File(void) {
  /* Function description:
  This function is designed to setup a simulation using the Setup.txt file.

  the setup file is case-insensitive. Thus, the reader is case-insensitive. */

  // Open setup.txt
  std::ifstream File;
  File.open("./IO/Setup.txt");
  if(File.is_open() == false) {
    throw Cant_Open_File("Cant Open File exception: Thrown by Simulation::Load_Setup_File\n"
                         "For some reason, ./IO/Setup.txt could not be opened\n");
  } // if(File.is_open() == false) {

  // Set up buffers
  std::string strBuf;

  // Now read simulation parameters

  //////////////////////////////////////////////////////////////////////////////
  // IO
  strBuf = IO::read_line_after(File, "Load Simulation from Save:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Load_Simulation_From_Save = true; }
  else {                                                        Simulation::Load_Simulation_From_Save = false; }

  strBuf = IO::read_line_after(File, "Save Data To File:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Save_Simulation_To_File = true; }
  else {                                                        Simulation::Save_Simulation_To_File = false; }

  strBuf = IO::read_line_after(File, "Print Particle Forces:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Print_Particle_Forces = true; }
  else {                                                        Simulation::Print_Particle_Forces = false; }

  strBuf = IO::read_line_after(File, "Print Body Forces:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Print_Body_Forces = true; }
  else {                                                        Simulation::Print_Body_Forces = false; }

  strBuf = IO::read_line_after(File, "Print Body Torques:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Print_Body_Torques = true; }
  else {                                                        Simulation::Print_Body_Torques = false; }

  strBuf = IO::read_line_after(File, "Print Box Boundary Forces:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::Print_Box_Boundary_Forces = true; }
  else {                                                        Simulation::Print_Box_Boundary_Forces = false; }

  strBuf = IO::read_line_after(File, "Time Steps Between Prints:", false);
  sscanf(strBuf.c_str(), " %u \n", &Simulation::TimeSteps_Between_Prints);


  #if defined(SIMULATION_SETUP_MONITOR)
    printf(  "Read Load_Simulation_From_Save as:      ");
    if(Simulation::Load_Simulation_From_Save == true) {    printf("true\n");  }
    else {                                                 printf("false\n"); }

    printf(  "Read Save_Simulation_To_File as:        ");
    if(Simulation::Save_Simulation_To_File == true) {      printf("true\n");  }
    else {                                                 printf("false\n"); }

    printf(  "Read Print_Particle_Forces as:          ");
    if(Simulation::Print_Particle_Forces == true) {        printf("true\n");  }
    else {                                                 printf("false\n"); }

    printf(  "Read Print_Body_Forces as:              ");
    if(Simulation::Print_Body_Forces == true) {            printf("true\n");  }
    else {                                                 printf("flase\n"); }

    printf(  "Read Print_Body_Torques as:             ");
    if(Simulation::Print_Body_Torques == true) {           printf("true\n");  }
    else {                                                 printf("flase\n"); }

    printf(  "Read Print_Box_Boundary_Forces as:      ");
    if(Simulation::Print_Box_Boundary_Forces == true) {    printf("true\n");  }
    else {                                                 printf("flase\n"); }

    printf(  "Read Time_Steps_Between_Prints as:      ");
    if(Simulation::TimeSteps_Between_Prints == true) {     printf("true\n");  }
    else {                                                 printf("false\n"); }
  #endif


  //////////////////////////////////////////////////////////////////////////////
  // Time Step
  strBuf = IO::read_line_after(File, "dt [s]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::dt);

  strBuf = IO::read_line_after(File, "Num Time Steps:", false);
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Time_Steps);



  //////////////////////////////////////////////////////////////////////////////
  // Contact
  strBuf = IO::read_line_after(File, "Contact Distance [mm]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Contact_Distance);

  strBuf = IO::read_line_after(File, "Friction Coefficient:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Friction_Coefficient);

  #if defined(SIMULATION_SETUP_MONITOR)
    printf(  "Read TimeSteps_Between_Prints as:       %u\n",  Simulation::TimeSteps_Between_Prints);
    printf(  "Read dt as:                             %lf (s)\n", Simulation::dt);
    printf(  "Read Num_Time_Steps as:                 %u\n",  Simulation::Num_Time_Steps);
    printf(  "Read Contact_Distance as:               %lf (mm)\n", Simulation::Contact_Distance);
    printf(  "Read Friction_Coefficient as:           %lf\n", Simulation::Friction_Coefficient);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // If "Load_Simulation_From_Save" is true, then stop reading the setup file.
  // Otherwise, read in the number of bodies
  if(Simulation::Load_Simulation_From_Save == true) {
    File.close();
    return nullptr;
  } // if(Simulation::Load_Simulation_From_Save == true) {




  //////////////////////////////////////////////////////////////////////////////
  // Number of Bodies
  strBuf = IO::read_line_after(File, "Number of Bodies:", false);
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Bodies);

  #if defined(SIMULATION_SETUP_MONITOR)
    printf(  "Read Number of Bodies as:               %u\n", Simulation::Num_Bodies);
  #endif

  // We're done reading in the Simulation parameters. We can close the setup file.
  File.close();

  // Allocate the bodies
  Body* Bodies = new Body[Simulation::Num_Bodies];

  // Now allocate Body-specific members.
  Simulation::From_FEB_File            = new bool[Num_Bodies];
  Simulation::General_BCs              = new Array<General_Boundary_Condition>[Num_Bodies];
  Simulation::Position_Offset          = new Vector[Num_Bodies];
  Simulation::Initial_Velocity         = new Vector[Num_Bodies];

  // Load the bodies from the setup file.
  for(unsigned i = 0; i < Simulation::Num_Bodies; i++) { Load_Body_From_Setup_File(Bodies[i], i); }

  // All done!
  return Bodies;
} // Body* Simulation::Load_Setup_File(void) {





void Simulation::Load_Body_From_Setup_File(Body & Body_In, const unsigned i) {
  /* Function description:
  This function reads in the ith body from Setup.txt.

  Load_Setup_File is the only thing that should call this function. */

  // Open setup.txt
  std::ifstream File;
  File.open("./IO/Setup.txt");
  if(File.is_open() == false) {
    throw Cant_Open_File("Cant Open File exception: Thrown by Simulation::Load_Body_From_Setup_File\n"
                         "For some reason, ./IO/Setup.txt could not be opened\n");
  } // if(File.is_open() == false) {

  // Set up buffers
  std::string strBuf;
  unsigned uBuf1, uBuf2, uBuf3;
  double lfBuf;
  char Buf[256];
  Simulation::Box_BCs Box_Boundary_Conditions;

  // First, go to the line of specified body.
  strBuf = "# Body ";
  strBuf += std::to_string(i);
  IO::read_line_after(File, strBuf.c_str(), false);



  //////////////////////////////////////////////////////////////////////////////
  // General

  strBuf = IO::read_line_after(File, "Name:", false);
  sscanf(strBuf.c_str(), " %s \n", Buf);
  Body_In.Set_Name(std::string{Buf});

  strBuf = IO::read_line_after(File, "From FEB File:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Simulation::From_FEB_File[i] = true; }
  else {                                                        Simulation::From_FEB_File[i] = false; }

  bool Is_Box;
  strBuf = IO::read_line_after(File, "Is Box:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Is_Box = true; }
  else {                                                        Is_Box = false; }

  // Box specific junk
  if(Is_Box == true) {
    strBuf = IO::read_line_after(File, "Dimensions [# Particles]:", false);
    sscanf(strBuf.c_str(), " {%u, %u, %u} \n", &uBuf1, &uBuf2, &uBuf3);
    Body_In.Set_Box_Dimensions(uBuf1, uBuf2, uBuf3);

    strBuf = IO::read_line_after(File, "x plus BC [mm/s]:", false);
    Box_Boundary_Conditions.x_plus_BC = Parse_BC(strBuf);

    strBuf = IO::read_line_after(File, "x minus BC [mm/s]:", false);
    Box_Boundary_Conditions.x_minus_BC = Parse_BC(strBuf);

    strBuf = IO::read_line_after(File, "y plus BC [mm/s]:", false);
    Box_Boundary_Conditions.y_plus_BC = Parse_BC(strBuf);

    strBuf = IO::read_line_after(File, "y minus BC [mm/s]:", false);
    Box_Boundary_Conditions.y_minus_BC = Parse_BC(strBuf);

    strBuf = IO::read_line_after(File, "z plus BC [mm/s]:", false);
    Box_Boundary_Conditions.z_plus_BC = Parse_BC(strBuf);

    strBuf = IO::read_line_after(File, "z minus BC [mm/s]:", false);
    Box_Boundary_Conditions.z_minus_BC = Parse_BC(strBuf);

    // Now set the BC
    Simulation::Set_Box_BCs(Body_In, Box_Boundary_Conditions);
  } // if(Is_Box[i] == true) {

  strBuf = IO::read_line_after(File, "Is Fixed In Place:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Body_In.Set_Is_Fixed(true); }
  else {                                                        Body_In.Set_Is_Fixed(false);  }

  strBuf = IO::read_line_after(File, "Is Damageable:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Body_In.Set_Is_Damageable(true); }
  else {                                                        Body_In.Set_Is_Damageable(false); }

  // If damageable then set damage rate.
  double Stretch_Critical_Mean, Stretch_Critical_SD;
  if(Body_In.Get_Is_Damageable() == true) {
    strBuf = IO::read_line_after(File, "Critical Stretch Mean:", false);
    sscanf(strBuf.c_str(), " %lf \n", &Stretch_Critical_Mean);

    strBuf = IO::read_line_after(File, "Critical Stretch SD:", false);
    sscanf(strBuf.c_str(), " %lf \n", &Stretch_Critical_SD);

    Body_In.Set_Particles_Critical_Stretch(Stretch_Critical_Mean, Stretch_Critical_SD);

    strBuf = IO::read_line_after(File, "Tau (Damage rate):", false);
    sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
    Body_In.Set_Tau(lfBuf);
  } // if(Body_In.Get_Is_Damageable() == true) {



  //////////////////////////////////////////////////////////////////////////////
  // Material

  Materials::Material Mat;


  // Read in the material parameters
  strBuf = IO::read_line_after(File, "Lame Parameter [Mpa]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Mat.Lame);

  strBuf = IO::read_line_after(File, "mu0 (Shear Modulus) [Mpa]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Mat.mu0);

  strBuf = IO::read_line_after(File, "Density [g/mm^3]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &Mat.density);


  // Now set the Body Material.
  Body_In.Set_Material(Mat);



  //////////////////////////////////////////////////////////////////////////////
  // Initial Conditions

  strBuf = IO::read_line_after(File, "Offset [mm]:", false);
  sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n", &Position_Offset[i][0], &Position_Offset[i][1], &Position_Offset[i][2]);

  strBuf = IO::read_line_after(File, "Initial Velocity [mm/s]:", false);
  sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n", &Initial_Velocity[i][0], &Initial_Velocity[i][1], &Initial_Velocity[i][2]);



  //////////////////////////////////////////////////////////////////////////////
  // General Boundary Conditions

  // First, read in number of BCs
  unsigned Number_General_BCs;
  strBuf = IO::read_line_after(File, "Number of General Boundary Conditions:", false);
  sscanf(strBuf.c_str(), " %u \n", &Number_General_BCs);

  // Check if there are any General BCs for this body
  if(Number_General_BCs != 0) {
    // If so, then make the BCs array
    General_BCs[i].Set_Length(Number_General_BCs);

    // Now read in each Boundary condition for this body
    for(unsigned j = 0; j < Number_General_BCs; j++) {
      strBuf = IO::read_line_after(File, "Condition Plane Normal Vector:", false);
      sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n",
            &General_BCs[i][j].Condition_Plane_Normal_Vector[0],
            &General_BCs[i][j].Condition_Plane_Normal_Vector[1],
            &General_BCs[i][j].Condition_Plane_Normal_Vector[2]);

      strBuf = IO::read_line_after(File, "Condition Plane Distance [mm]:", false);
      sscanf(strBuf.c_str(), " %lf \n", &General_BCs[i][j].Condition_Plane_Distance);

      /* Read in the condition inequality. Note that we must check for >= and <=
      BEFORE we check < and > (since the last two are contained in the first two) */
      strBuf = IO::read_line_after(File, "Condition Inequality:", false);
           if(IO::String_Ops::Contains(strBuf.c_str(), "<=")) { General_BCs[i][j].Condition_Inequality = Simulation::Inequality::LE; }
      else if(IO::String_Ops::Contains(strBuf.c_str(), ">=")) { General_BCs[i][j].Condition_Inequality = Simulation::Inequality::GE; }
      else if(IO::String_Ops::Contains(strBuf.c_str(), "==")) { General_BCs[i][j].Condition_Inequality = Simulation::Inequality::E;  }
      else if(IO::String_Ops::Contains(strBuf.c_str(), "<" )) { General_BCs[i][j].Condition_Inequality = Simulation::Inequality::L;  }
      else if(IO::String_Ops::Contains(strBuf.c_str(), ">" )) { General_BCs[i][j].Condition_Inequality = Simulation::Inequality::G;  }

      strBuf = IO::read_line_after(File, "Effect Vector [mm/s]:", false);
      sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n",
             &General_BCs[i][j].Effect_Vector[0],
             &General_BCs[i][j].Effect_Vector[1],
             &General_BCs[i][j].Effect_Vector[2]);
    } // for(unsigned j = 0; j < Number_General_BCs; j++) {
  } // if(Number_General_BCs != 0) {



  //////////////////////////////////////////////////////////////////////////////
  // Other

  strBuf = IO::read_line_after(File, "Enable Gravity:", false);
  if(IO::String_Ops::Contains(strBuf.c_str(), "true", false)) { Body_In.Set_Gravity_Enabled(true); }
  else {                                                        Body_In.Set_Gravity_Enabled(false); }

  strBuf = IO::read_line_after(File, "Mu (Viscosity) [Mpa*s]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_mu(lfBuf);

  strBuf = IO::read_line_after(File, "Alpha (Hg Control parameter):", false);
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_alpha(lfBuf);

  strBuf = IO::read_line_after(File, "Inter Particle Spacing (IPS) [mm]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_Inter_Particle_Spacing(lfBuf);

  strBuf = IO::read_line_after(File, "Support Radius [mm]:", false);
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_Support_Radius(lfBuf);

  #if defined(SIMULATION_SETUP_MONITOR)
    printf("\n\nRead Name[%3u] as:                      %s\n", i, Body_In.Get_Name().c_str());
    printf(    "Read From_FEB_File[%3u] as:             ", i);
    if(Simulation::From_FEB_File[i] == true) {   printf("true\n");  }
    else {                                       printf("false\n"); }

    printf("Read Is_Box[%3u] as:                    ", i);
    if(Body_In.Get_Is_Box() == true) {           printf("true\n");  }
    else {                                       printf("false\n"); }

    if(Is_Box == true) {
      printf("Read Dimensions[%3u] as:                {%u, %u, %u}\n", i, Body_In.Get_X_SIDE_LENGTH(), Body_In.Get_Y_SIDE_LENGTH(), Body_In.Get_Z_SIDE_LENGTH());
      printf("Read Body[%3u].x_plus_BC as:            ", i); Print_BC(Box_Boundary_Conditions.x_plus_BC);
      printf("Read Body[%3u].x_minus_BC as:           ", i); Print_BC(Box_Boundary_Conditions.x_minus_BC);
      printf("Read Body[%3u].y_plus_BC as:            ", i); Print_BC(Box_Boundary_Conditions.y_plus_BC);
      printf("Read Body[%3u].y_minus_BC as:           ", i); Print_BC(Box_Boundary_Conditions.y_minus_BC);
      printf("Read Body[%3u].z_plus_BC as:            ", i); Print_BC(Box_Boundary_Conditions.z_plus_BC);
      printf("Read Body[%3u].z_minus_BC as:           ", i); Print_BC(Box_Boundary_Conditions.z_minus_BC);
    } // if(Is_Box == true) {
    printf("Read Body[%3u].Is_Fixed as:             ", i);
    if(Body_In.Get_Is_Fixed() == true) {         printf("true\n");  }
    else {                                       printf("false\n"); }

    printf("Read Body[%3u].Is_Damageable as:        ", i);
    if(Body_In.Get_Is_Damageable() == true) {    printf("true\n");  }
    else {                                       printf("false\n"); }

    if(Body_In.Get_Is_Damageable() == true) {
      printf("Read Stretch Critical Mean as:          %lf\n", Stretch_Critical_Mean);
      printf("Read Stretch Critical SD as:            %lf\n", Stretch_Critical_SD);
      printf("Read Body[%3u].Tau as:                  %lf\n", i, Body_In.Get_Tau());
    } // if(Body_In.Get_Is_Damageable() == true) {

    printf("\nRead Body[%3u].Material.Lame as:        %lf (Mpa)\n", i, Mat.Lame);
    printf(  "Read Body[%3u].Material.mu0 as:         %lf (Mpa)\n", i, Mat.mu0);
    printf(  "Read Body[%3u].Material.density as:     %lf (g/mm^3)\n", i, Mat.density);

    printf("\nRead Offset[%3u] as:                    {%lf, %lf, %lf} (mm)\n",   i, Position_Offset[i][0], Position_Offset[i][1], Position_Offset[i][2]);
    printf(  "Read Initial_Velocity[%3u] as:          {%lf, %lf, %lf} (mm/s)\n", i, Initial_Velocity[i][0], Initial_Velocity[i][1], Initial_Velocity[i][2]);

    printf("\nRead Number_General_BCs[%3u] as:        %u\n",  i, Number_General_BCs);
    for(unsigned j = 0; j < Number_General_BCs; j++) {
      printf("Read Body[%3u]'s Plane_Vector %3u as:   {%lf, %lf, %lf}\n", i, j,
              General_BCs[i][j].Condition_Plane_Normal_Vector[0],
              General_BCs[i][j].Condition_Plane_Normal_Vector[1],
              General_BCs[i][j].Condition_Plane_Normal_Vector[2]);

      printf("Read Body[%3u]'s Plane_Distance %3u as: %lf\n", i, j, General_BCs[i][j].Condition_Plane_Distance);
      printf("Read Body[%3u]'s Inequality %3u as:     ", i, j);
      switch(General_BCs[i][j].Condition_Inequality) {
        case (Inequality::L):  printf("<\n" ); break;
        case (Inequality::LE): printf("<=\n"); break;
        case (Inequality::E):  printf("==\n"); break;
        case (Inequality::GE): printf(">=\n"); break;
        case (Inequality::G):  printf(">\n" ); break;
      } // switch(General_BCs[i][j].Condition_Inequality) {
      printf("Read Body[%3u]'s Effect_Vector %3u as:  {%lf, %lf, %lf}\n", i, j,
              General_BCs[i][j].Effect_Vector[0],
              General_BCs[i][j].Effect_Vector[1],
              General_BCs[i][j].Effect_Vector[2]);
    } // for(unsigned j = 0; j < Number_General_BCs; j++) {

    printf("\nRead Bodu[%3u].Gravity_Enabled          ", i);
    if(Body_In.Get_Gravity_Enabled() == true)    { printf("true\n");  }
    else                                         { printf("false\n"); }

    printf(  "Read Body[%3u].mu as:                   %lf (Mpa*s)\n", i, Body_In.Get_mu());
    printf(  "Read Body[%3u].alpha as:                %lf\n",         i, Body_In.Get_alpha());
    printf(  "Read Body[%3u].IPS as:                  %lf (mm)\n",    i, Body_In.Get_Inter_Particle_Spacing());
    printf(  "Read Body[%3u].Support_Radius as:       %lf (mm)\n",    i, Body_In.Get_Support_Radius());
  #endif


  // All done. Close file
  File.close();
} // void Simulation::Load_Body_From_Setup_File(Body & Body_In, const unsigned i) {





Vector Simulation::Parse_BC(std::string & Line) {
  /* Function description:
  This function is designed to take a line of the form " {BCx, BCy, BCz} and
  return the vector <BCx, BCy, BCz>, where BCx/y/z is either a double or "free".

  Load_Body_From_Setup_File is the only thing that should call this function. */

  // Setup buffers
  char Buf[3][32];

  // Read in the 3 components of the BC, store them in the buffers.
  sscanf(Line.c_str(), " {%s %s %s} \n", Buf[0], Buf[1], Buf[2]);

  // Ensure that Buf1/2/3 are null-terminated
  Buf[0][31] = Buf[1][31] = Buf[2][31] = '\0';

  // Set up Vector that we'll return
  Vector Parsed_BC{};

  /* For each buffer, check if it contains "Free". If so, then set the
  corresponding component of Parsed_BC to Simulation::FREE. Otherwise,
  set the corresponding component to the value that we read in. */
  for(unsigned i = 0; i < 3; i++) {
    if(IO::String_Ops::Contains(Buf[i], "Free", false)) { Parsed_BC[i] = Simulation::FREE; }
    else { sscanf(Buf[i], "%lf", &Parsed_BC[i]); }
  } // for(unsigned i = 0; i < 3; i++) {

  return Parsed_BC;
} // void Simulation::Parse_BC(std::string & Line) {



void Simulation::Print_BC(Vector & BC_In) {
  /* Function description:
  This function is for testing only. It exists only to eliminate code repetition */
  printf("{");
  for(unsigned j = 0; j < 3; j++) {
    if(BC_In[j] == Simulation::FREE) { printf("Free"); }
    else { printf("%lf", BC_In[j]); }

    if(j == 2) { printf("}\n"); }
    else { printf(", "); }
  } // for(unsigned j = 0; j < 3; j++) {
} // void Simulation::Print_BC(Vector & BC_In) {
