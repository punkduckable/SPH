#include "Simulation.h"
#include "IO/IO_Ops.h"
#include "Errors.h"
#include "Body/Body.h"
#include <fstream>
#include <string>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
// Define external simulation.h variables

// Body properties
namespace Simulation {
  bool Load_Simulation_From_Save           = true;
  bool Save_Simulation_To_File             = false;
  bool Print_Particle_Forces               = false;
  bool Print_Net_External_Forces           = true;
  unsigned TimeSteps_Between_Prints        = 1000;


  // TimeStep paramters
  double dt                                = .0000001;        // Time step        : s
  unsigned Num_Time_Steps                  = 10000;           // Number of time steps

  // Contact
  double Contact_Distance = 1;                  // Distance at which bodies begin contacting one another.   : mm
  double Friction_Coefficient = .1;                                            //        : unitless

  unsigned Num_Bodies = 0;                       // Number of bodies in simulation

  std::string * Names = nullptr;                 // The names of each body (name must match File name if reading from FEB file)
  bool * Is_Box = nullptr;                       // Which bodies are Boxs
  Box_Properties * Box_Parameters = nullptr;     // Specifies the dimensions, and BCs of the box bodies.
  bool * Is_Fixed = nullptr;                     // Which bodies are fixed in place (can be from FEB file or Box)
  bool * Is_Damagable = nullptr;                 // Which bodies can be damaged
  bool * From_FEB_File = nullptr;                // Which bodies will be read from file
  unsigned * Time_Steps_Per_Update=nullptr;      // How many time steps pass between updating this Body's P-K tensor
  double * IPS = nullptr;                        // Inter particle spacing in mm.
  Vector * Position_Offset = nullptr;            // Position offset for particles in body
  Vector * Initial_Velocity = nullptr;           // Initial velocity condition
  Materials::Material * Simulation_Materials = nullptr;    // Each bodies material
} // namespace Simulation {



Body* Simulation::Load_Setup_File(void) {
  /* Function description:
  This function is designed to setup a simulation using the Setup.txt file. */

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
  strBuf = IO::read_line_after(File, "Load Simulation from Save:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) { Load_Simulation_From_Save = true; }
  else {                                                 Load_Simulation_From_Save = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read Load_Simulation_From_Save as: %s\n", strBuf.c_str());
  #endif

  strBuf = IO::read_line_after(File, "Save Data To File:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) { Save_Simulation_To_File = true; }
  else {                                                 Save_Simulation_To_File = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read Save_Simulation_To_File as: %s\n", strBuf.c_str());
  #endif

  strBuf = IO::read_line_after(File, "Print Particle Forces:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) { Print_Particle_Forces = true; }
  else {                                                 Print_Particle_Forces = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read Print_Particle_Forces as: %s\n", strBuf.c_str());
  #endif

  strBuf = IO::read_line_after(File, "Print Net External Forces:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) { Print_Net_External_Forces = true; }
  else {                                                 Print_Net_External_Forces = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read Print_Net_External_Forces as: %s\n", strBuf.c_str());
  #endif

  strBuf = IO::read_line_after(File, "Time Steps Between Prints:");
  sscanf(strBuf.c_str(), " %u \n", &Simulation::TimeSteps_Between_Prints);
  #if defined(SIMULATION_MONITOR)
    printf("Read TimeSteps_Between_Prints as:       %u\n", Simulation::TimeSteps_Between_Prints);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Time Step
  strBuf = IO::read_line_after(File, "dt:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::dt);
  #if defined(SIMULATION_MONITOR)
    printf("Read dt as:                             %lf\n", Simulation::dt);
  #endif

  strBuf = IO::read_line_after(File, "Num Time Steps:");
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Time_Steps);
  #if defined(SIMULATION_MONITOR)
    printf("Read Num_Time_Steps as:                 %u\n", Simulation::Num_Time_Steps);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Contact
  strBuf = IO::read_line_after(File, "Contact Distance:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Contact_Distance);
  #if defined(SIMULATION_MONITOR)
    printf("Read Contact_Distance as:               %lf\n", Simulation::Contact_Distance);
  #endif

  strBuf = IO::read_line_after(File, "Friction Coefficient:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Friction_Coefficient);
  #if defined(SIMULATION_MONITOR)
    printf("Read Friction_Coefficient as:           %lf\n", Simulation::Friction_Coefficient);
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
  strBuf = IO::read_line_after(File, "Number of Bodies:");
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Bodies);
  #if defined(SIMULATION_MONITOR)
    printf("Read Load_Simulation_From_Save as:      %u\n", Simulation::Num_Bodies);
  #endif

  // We're done reading in the Simulation parameters. We can close the setup file.
  File.close();

  // Allocate the bodies
  Body * Bodies = new Body[Simulation::Num_Bodies];

  // Now allocate Body-specific members.
  From_FEB_File = new bool[Num_Bodies];
  Box_Parameters = new Box_Properties[Num_Bodies];
  Position_Offset = new Vector[Num_Bodies];
  Initial_Velocity = new Vector[Num_Bodies];


  // Load the bodies from the setup file.
  for(unsigned i = 0; i < Simulation::Num_Bodies; i++) { Load_Body_From_Setup_File(Bodies[i], i); }

  return Bodies;
} // Body* Simulation::Load_Setup_File(void) {





void Simulation::Load_Body_From_Setup_File(Body & Body_In, const unsigned i) {
  /* Function description:
  This function reads in the ith body from Setup.txt. */
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

  // First, go to the line of specified body.
  strBuf = "# Body ";
  strBuf += std::to_string(i);
  IO::read_line_after(File, strBuf.c_str());



  //////////////////////////////////////////////////////////////////////////////
  // General

  strBuf = IO::read_line_after(File, "Name:");
  sscanf(strBuf.c_str(), " %s \n", Buf);
  Body_In.Set_Name(std::string{Buf});
  #if defined(SIMULATION_MONITOR)
    printf("\nRead Name[%3u] as:                      %s\n", i, Body_In.Get_Name().c_str());
  #endif

  strBuf = IO::read_line_after(File, "From FEB File:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) {      Simulation::From_FEB_File[i] = true; }
  else {                                                      Simulation::From_FEB_File[i] = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read From_FEB_File[%3u] as:             %u\n", i, Simulation::From_FEB_File[i]);
  #endif

  bool Is_Box;
  strBuf = IO::read_line_after(File, "Is Box:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) {      Is_Box = true; }
  else {                                                      Is_Box = false; }
  #if defined(SIMULATION_MONITOR)
    printf("Read Is_Box[%3u] as:                    %u\n", i, Is_Box);
  #endif

  // Box specific junk
  if(Is_Box == true) {
    strBuf = IO::read_line_after(File, "Dimensions:");
    sscanf(strBuf.c_str(), " {%u, %u, %u} \n", &uBuf1, &uBuf2, &uBuf3);
    Body_In.Set_Box_Dimensions(uBuf1, uBuf2, uBuf3);
    #if defined(SIMULATION_MONITOR)
      printf("Read Dimensions[%3u] as:                {%u, %u, %u}\n", i, uBuf1, uBuf2, uBuf3);
    #endif


    // BCs
  } // if(Is_Box[i] == true) {

  strBuf = IO::read_line_after(File, "Is Fixed In Place:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) {      Body_In.Set_Is_Fixed(true); }
  else {                                                      Body_In.Set_Is_Fixed(false);  }
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Is_Fixed as:             %u\n", i, Body_In.Get_Is_Fixed());
  #endif

  strBuf = IO::read_line_after(File, "Is Damageable:");
  if(IO::String_Ops::Contains(strBuf.c_str(), "true")) {      Body_In.Set_Is_Damageable(true); }
  else {                                                      Body_In.Set_Is_Damageable(false); }
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Is_Damageable as:        %u\n", i, Body_In.Get_Is_Damageable());
  #endif

  // If damageable then set damage rate.
  if(Body_In.Get_Is_Damageable() == true) {
    strBuf = IO::read_line_after(File, "Tau (Damage rate):");
    sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
    Body_In.Set_Tau(lfBuf);
    #if defined(SIMULATION_MONITOR)
      printf("Read Body[%3u].Tau as:                  %lf\n", i, Body_In.Get_Tau());
    #endif
  } // if(Body_In.Get_Is_Damageable() == true) {


  /*
      x_plus_BC:                                   {0, Free, Free}
      x_minus_BC:                                  {0, Free, Free}
      y_plus_BC:                                   {Free, 0, Free}
      y_minus_BC:                                  {Free, 0, Free}
      z_plus_BC:                                   {Free, Free, 0}
      z_minus_BC:                                  {Free, Free, 0}
    */



  //////////////////////////////////////////////////////////////////////////////
  // Material

  Materials::Material Mat;

  strBuf = IO::read_line_after(File, "Lame Parameter (Mpa):");
  sscanf(strBuf.c_str(), " %lf \n", &Mat.Lame);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Material.Lame as:        %lf\n", i, Mat.Lame);
  #endif

  strBuf = IO::read_line_after(File, "mu0 (Shear Modulus) (Mpa):");
  sscanf(strBuf.c_str(), " %lf \n", &Mat.mu0);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Material.mu0 as:         %lf\n", i, Mat.mu0);
  #endif

  strBuf = IO::read_line_after(File, "Density (g/mm^3):");
  sscanf(strBuf.c_str(), " %lf \n", &Mat.density);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Material.density as:     %lf\n", i, Mat.density);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Initial Conditions

  strBuf = IO::read_line_after(File, "Offset:");
  sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n", &Position_Offset[i][0], &Position_Offset[i][1], &Position_Offset[i][2]);
  #if defined(SIMULATION_MONITOR)
    printf("Read Offset[%3u] as:                    {%lf, %lf, %lf}\n", i, Position_Offset[i][0], Position_Offset[i][1], Position_Offset[i][2]);
  #endif

  strBuf = IO::read_line_after(File, "Initial Velocity:");
  sscanf(strBuf.c_str(), " {%lf, %lf, %lf} \n", &Initial_Velocity[i][0], &Initial_Velocity[i][1], &Initial_Velocity[i][2]);
  #if defined(SIMULATION_MONITOR)
    printf("Read Initial_Velocity[%3u] as:          {%lf, %lf, %lf}\n", i, Initial_Velocity[i][0], Initial_Velocity[i][1], Initial_Velocity[i][2]);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Other

  strBuf = IO::read_line_after(File, "Mu (Viscosity):");
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_mu(lfBuf);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].mu as:                   %lf\n", i, Body_In.Get_mu());
  #endif

  strBuf = IO::read_line_after(File, "Alpha (Hg Control parameter):");
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_alpha(lfBuf);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].alpha as:                %lf\n", i, Body_In.Get_alpha());
  #endif

  strBuf = IO::read_line_after(File, "Support Radius:");
  sscanf(strBuf.c_str(), " %u \n", &uBuf1);
  Body_In.Set_Support_Radius(uBuf1);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Support_Radius as:       %u\n", i, Body_In.Get_Support_Radius());
  #endif

  strBuf = IO::read_line_after(File, "Inter Particle Spacing:");
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_Inter_Particle_Spacing(lfBuf);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].IPS as:                  %lf\n", i, Body_In.Get_Inter_Particle_Spacing());
  #endif

  strBuf = IO::read_line_after(File, "Time Steps Per Update:");
  sscanf(strBuf.c_str(), " %u \n", &uBuf1);
  Body_In.Set_Time_Steps_Per_Update(uBuf1);
  #if defined(SIMULATION_MONITOR)
    printf("Read Body[%3u].Time_Steps/Update as:    %u\n", i, Body_In.Get_Time_Steps_Per_Update());
  #endif

  // All done. Close file
  File.close();
} // void Simulation::Load_Body_From_Setup_File(Body & Body_In, const unsigned i) {
