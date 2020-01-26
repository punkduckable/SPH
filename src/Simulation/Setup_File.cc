#include "Simulation.h"
#include "IO/IO_Ops.h"
#include "Errors.h"
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



void Simulation::Load_Setup_File(void) {
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
    printf("Read TimeSteps_Between_Prints as:    %u\n", Simulation::TimeSteps_Between_Prints);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Time Step
  strBuf = IO::read_line_after(File, "dt:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::dt);
  #if defined(SIMULATION_MONITOR)
    printf("Read dt as:                          %lf\n", Simulation::dt);
  #endif

  strBuf = IO::read_line_after(File, "Num Time Steps:");
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Time_Steps);
  #if defined(SIMULATION_MONITOR)
    printf("Read Num_Time_Steps as:              %u\n", Simulation::Num_Time_Steps);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // Contact
  strBuf = IO::read_line_after(File, "Contact Distance:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Contact_Distance);
  #if defined(SIMULATION_MONITOR)
    printf("Read Contact_Distance as:            %lf\n", Simulation::Contact_Distance);
  #endif

  strBuf = IO::read_line_after(File, "Friction Coefficient:");
  sscanf(strBuf.c_str(), " %lf \n", &Simulation::Friction_Coefficient);
  #if defined(SIMULATION_MONITOR)
    printf("Read Friction_Coefficient as:        %lf\n", Simulation::Friction_Coefficient);
  #endif



  //////////////////////////////////////////////////////////////////////////////
  // If "Load_Simulation_From_Save" is true, then stop reading the setup file.
  // Otherwise, read in the number of bodies
  if(Simulation::Load_Simulation_From_Save == true) {
    File.close();
    return;
  } // if(Simulation::Load_Simulation_From_Save == true) {




  //////////////////////////////////////////////////////////////////////////////
  // Number of Bodies
  strBuf = IO::read_line_after(File, "Number of Bodies:");
  sscanf(strBuf.c_str(), " %u \n", &Simulation::Num_Bodies);
  #if defined(SIMULATION_MONITOR)
    printf("Read Load_Simulation_From_Save as:   %u\n", Simulation::Num_Bodies);
  #endif

  // We're done reading in the Simulation parameters. We can close the setup file.
  File.close();

  // Now allocate Body-specific members.
  Names = new std::string[Num_Bodies];
  From_FEB_File = new bool[Num_Bodies];
  Is_Box = new bool[Num_Bodies];
  Box_Parameters = new Box_Properties[Num_Bodies];
  Is_Fixed = new bool[Num_Bodies];
  Is_Damagable = new bool[Num_Bodies];
  Time_Steps_Per_Update = new unsigned[Num_Bodies];
  IPS = new double[Num_Bodies];
  Position_Offset = new Vector[Num_Bodies];
  Initial_Velocity = new Vector[Num_Bodies];


  // Load the bodies from the setup file.
  for(unsigned i = 0; i < Simulation::Num_Bodies; i++) { Load_Body_From_Setup_File(i); }
} // void Simulation::Load_Setup_File(void) {


void Simulation::Load_Body_From_Setup_File(const unsigned i) {
  /* Function description:
  This function reads in the ith body from Setup.txt. */
  // Open setup.txt
  std::ifstream File;
  File.open("./IO/Setup.txt");
  if(File.is_open() == false) {
    throw Cant_Open_File("Cant Open File exception: Thrown by Simulation::Load_Body_From_Setup_File\n"
                         "For some reason, ./IO/Setup.txt could not be opened\n");
  } // if(File.is_open() == false) {

  // Set up buffer
  std::string strBuf;

  // Close file
  File.close();
} // void Simulation::Load_Body_From_Setup_File(const unsigned i) {
