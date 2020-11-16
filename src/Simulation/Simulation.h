#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

//#define SIMULATION_DEBUG
#define SIMULATION_SETUP_MONITOR

#include "Classes.h"
#include "Vector/Vector.h"
#include "Materials.h"
#include <string>
#include <ctime>

const double CLOCKS_PER_MS = CLOCKS_PER_SEC/1000.;

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  //////////////////////////////////////////////////////////////////////////////
  // Setup functions.
  // Defined in Simulation_Setup.cc

  void Setup(Body ** Bodies);


  //////////////////////////////////////////////////////////////////////////////
  // Functions that run simulations
  // Defined in Simulation.cc

  void Run(void);


  //////////////////////////////////////////////////////////////////////////////
  // Timing functions
  // Defined in Timing.cc

  #if defined(_OPENMP)
    #define TIME_TYPE double
  #else
    #define TIME_TYPE clock_t
  #endif
  TIME_TYPE Get_Time(void);
  TIME_TYPE Time_Since(TIME_TYPE time);


  //////////////////////////////////////////////////////////////////////////////
  // Setup file function
  // Defined in Setup_File.cc

  Body* Load_Setup_File(void);


  //////////////////////////////////////////////////////////////////////////////
  // Simulation parameters
  // Declared in Setup_File.cc (defined by Load_Setup_File.cc and Setup.txt)

  // IO paramaters
  extern bool Load_Simulation_From_Save;
  extern bool Save_Simulation_To_File;
  extern bool Print_Particle_Forces;
  extern bool Print_Body_Forces;
  extern unsigned TimeSteps_Between_Prints;

  // Time Step
  extern double dt;                              // Time step                   : s
  extern unsigned Num_Time_Steps;                // Number of time steps

  // Contact
  extern double Contact_Distance;                // Distance at which bodies begin contacting one another.   : mm
  extern double Friction_Coefficient;                                          // unitless

  // Number of bodies
  extern unsigned Num_Bodies;                    // Number of bodies in simulation


  //////////////////////////////////////////////////////////////////////////////
  // Simulation Setup Parameters
  // These variables are used to set up the simulation. They are defined in
  // Setup_File.cc and used in Simulation_Setup.cc

  /* To simplify implementation, I choose a random value to designate as "free".
  If a particular BC component has this value, then that BC component is treated
  as if it were "free" (had no BC applied to it). This means that it is impossible
  to set a box BC of exactly -293103918 mm/s */
  const double FREE = -293103918;

  enum class Inequality{L = -2, LE = -1, E = 0, GE = 1, G = 2};
  struct General_Boundary_Condition {
    /* This structure is used to define a General Boundary condition. BCs affect
    the velocity of particles in a body.
    Each BC consists of two parts: a condition and an effect.


    Condition: The condition is used to determine which particles (in a body)
    will have the BC applied to them. The condition consists of....
        Normal_Plane_Vector: Normal vector of a plane in 3d space

        Distance_To_Plane: The distance, in the direction of the Normal_Vector,
        from the origin to the plane. The plane consists of all points in 3d
        space such that
            <x,y,z> dot Condition_Normal_Vector = Condition_Distance_To_Plane

        Inequality: Defines how the BC is applied relative to the plane
        (LE = Less than or equal, GE = Greater than or Equal).

    To understand how this works, suppose that
        Condition_Plane_Normal_Vector = {1,2,3}
        Condition_Plane_Distance = 15
        Condition_Inequality = GE
    Then any particle P whose reference position P_X satisifies
        Dot_Product(P_X, {1,2,3}) / Magnitude({1,2,3}) >= 15
    will have the BC effect applied to it.


    Effect: This determines the effect that is applied to all particles (in the
    body) that satisify the condition. The Effect consists of the Effect_Vector,
    which determines how the BC affects the velocity of the particles that
    satisify the condition. If you want to apply a BC in a particular direction,
    then the corresponding component of the Effect_Vector should have that value.
    If you do not want to apply a BC in a particular driection, then the
    corresponding component of the Effect_Vector should be "FREE" (which is a
    constant that's defined above. This means that you can not set a component
    of the Effect_Vector to the value of the FREE constant. If you do, the code
    will treat that component as Free/prescribe no BC. */

    // Condition
    Vector                   Condition_Plane_Normal_Vector;
    double                   Condition_Plane_Distance;
    Inequality               Condition_Inequality;

    // Effect
    Vector                   Effect_Vector;
  }; // struct General_Boundary_Condition {

  extern bool * From_FEB_File;                             // Which bodies will be read from file
  extern Array<General_Boundary_Condition> * General_BCs;  // Specifies the general boundary conditions for each body.
  extern Vector * Position_Offset;                         // Position offset for particles in body
  extern Vector * Initial_Velocity;                        // Initial velocity condition for each body



  //////////////////////////////////////////////////////////////////////////////
  // Boundary Conditions
  // Defined in Boundary_Conditions.cc

  struct Box_BCs {
    Vector x_plus_BC;
    Vector x_minus_BC;
    Vector y_plus_BC;
    Vector y_minus_BC;
    Vector z_plus_BC;
    Vector z_minus_BC;
  }; // struct Box_BCs {

  void Set_General_BCs(Body & Body_In,                     // The body we're applying the BC to
                       Array<General_Boundary_Condition> & BCs_In);  // The BCs being applied
  void Set_Box_BCs(Body & Box,                             // Reference to the box body
                   Box_BCs & Boundary_Conditions);         // Box's parameters
  void Set_Box_Particle_BCs(Particle & P_In,               // Particle that we're applying the BC to
                            Vector BC);                    // The BC that's being applied
} // namespace Simulation {

#endif
