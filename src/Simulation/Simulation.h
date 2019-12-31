#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#include "Classes.h"
#include "Materials.h"
#include <string>

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  //////////////////////////////////////////////////////////////////////////////
  // Setup functions.
  // Defined in Simulation_Setup.cc
  void Setup_Cuboid(Body & Body_In, const unsigned m);
  void Setup_FEB_Body(Body & FEB_Body, const unsigned m);
  void Body_Needle_Set_Up(void);                           // Set up Body/Needle simulation
  void Set_Body_Members(Body & Body_In);                   // Set default body members



  //////////////////////////////////////////////////////////////////////////////
  // Run the simulation.
  // Defined in Simulation.cc
  void Run_Simulation(void);



  //////////////////////////////////////////////////////////////////////////////
  // Simulation parameters
  // (Note: const variables have internal linkage)

  // Simulation flags/properties. Defined
  const unsigned char Load_Data_From_File        = 0;
  const unsigned char Save_Data_To_File          = 0;
  const unsigned char Print_Forces               = 1;
  const unsigned char Print_Net_Force            = 0;
  const unsigned TimeSteps_Between_Prints        = 1000;

  // TimeStep paramters
  const double dt                                = 5e-9;            // Time step        : s
  const unsigned Num_Steps                       = 100000;          // Number of time steps

  // Contact parameter
  const double Contact_Distance = 1;             // Distance at which bodies begin contacting one another.   : mm

  // Friction coefficient.
  const double Friction_Coefficient = .1;                                      //        : unitless

  // Body properties (defined in Simulation_Setup.cc)
  extern unsigned Num_Arrays;                    // Number of bodies in simulation
  extern std::string * Names;                    // The names of each body (name must match File name if reading from FEB file)
  extern bool * Is_Cuboid;                       // Which bodies are cuboids
  extern bool * Is_Boundary;                     // Which bodies are boundaries (can be from FEB file or cuboid)
  extern bool * Is_Damagable;                    // Which bodies can be damaged
  extern bool * From_FEB_File;                   // Which bodies will be read from file
  extern unsigned * Steps_Per_Update;            // How many time steps pass between updating this Body's P-K tensor
  extern double * IPS;                           // Inter particle spacing in mm.
  extern Vector * Dimensions;                    // Dimensions of cuboids (only applicable for cuboids)
  extern Vector * Offset;                        // Poisition offset (only applicable for cuboids)
  extern Vector * Initial_Velocity;              // Initial velocity condition
  extern double * Num_Particles;                 // The number of particles in each body
  extern Materials::Material * Simulation_Materials;       // Each bodies material
} // namespace Simulation {

#endif
