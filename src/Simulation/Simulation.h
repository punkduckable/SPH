#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#include "Classes.h"
#include "Materials.h"
#include <string>

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  void Run_Simulation(void);
  void Setup_Cuboid(Body & Body, const unsigned m);
  void Setup_FEB_Body(Body & FEB_Body, const unsigned m);

  // Simulation flags/properties
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

  void Use_Arrays_From_Code(void);

  // Default Body members
  void Set_Body_Members(Body & Body_In);
} // namespace Simulation {

#endif
