#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

//#define SIMULATION_DEBUG
#define SIMULATION_SETUP_MONITOR

#include "Classes.h"
#include "Vector/Vector.h"
#include "Materials.h"
#include <string>

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  //////////////////////////////////////////////////////////////////////////////
  // Setup functions.
  // Defined in Simulation_Setup.cc

  void Setup_Box(Body & Body_In, const unsigned m);
  void Setup_FEB_Body(Body & FEB_Body, const unsigned m);
  void Bodies_Setup(void);                                 // Set up Body/Needle simulation
  void Set_Body_Members(Body & Body_In);                   // Set default body members
  void Setup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index);


  //////////////////////////////////////////////////////////////////////////////
  // Functions that run simulations
  // Defined in Simulation.cc

  void Run_Simulation(void);


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

  // Simulation flags/properties. Defined
  extern bool Load_Simulation_From_Save;
  extern bool Save_Simulation_To_File;
  extern bool Print_Particle_Forces;
  extern bool Print_Net_External_Forces;
  extern unsigned TimeSteps_Between_Prints;

  // TimeStep paramters
  extern double dt;                              // Time step                   : s
  extern unsigned Num_Time_Steps;                // Number of time steps

  // Contact
  extern double Contact_Distance;                // Distance at which bodies begin contacting one another.   : mm
  extern double Friction_Coefficient;                                          // unitless



  //////////////////////////////////////////////////////////////////////////////
  // Body properties (defined in Simulation_Setup.cc)

  extern unsigned Num_Bodies;                    // Number of bodies in simulation

  /* To simplify implementation, I choose a random value to designate as "free".
  If a particular BC component has this value, then that BC component is treated
  as if it were "free" (had no BC applied to it). This means that it is impossible
  to set a box BC of exactly -293103918 mm/s */
  const double FREE = -293103918;

  struct Box_Properties {
    Vector Dimensions;

    // BCs
    Vector x_plus_BC;
    Vector x_minus_BC;
    Vector y_plus_BC;
    Vector y_minus_BC;
    Vector z_plus_BC;
    Vector z_minus_BC;
  }; // struct Box_Properties {

  struct Box_BCs {
    Vector x_plus_BC;
    Vector x_minus_BC;
    Vector y_plus_BC;
    Vector y_minus_BC;
    Vector z_plus_BC;
    Vector z_minus_BC;
  }; // struct Box_BCs {

  extern std::string * Names;                    // The names of each body (name must match File name if reading from FEB file)
  extern bool * Is_Box;                          // Which bodies are Boxs
  extern bool * Is_Fixed;                        // Which bodies are fixed in place (can be from FEB file or Box)
  extern bool * Is_Damagable;                    // Which bodies can be damaged
  extern bool * From_FEB_File;                   // Which bodies will be read from file
  extern unsigned * Time_Steps_Per_Update;       // How many time steps pass between updating this Body's P-K tensor
  extern double * IPS;                           // Inter particle spacing in mm.
  extern Box_Properties * Box_Parameters;        // Specifies the dimensions, and BCs of the box bodies.
  extern Vector * Position_Offset;               // Position offset for particles in body
  extern Vector * Initial_Velocity;              // Initial velocity condition
  extern Materials::Material * Simulation_Materials;       // Each bodies material



  //////////////////////////////////////////////////////////////////////////////
  // Boundary Conditions
  // Defined in Boundary_Conditions.cc

  enum class Position_Type{Reference = 0, Spatial = 1};
  enum class Inequality{LE = -1, GE = 1};

  struct Boundary_Condition {
    /* This structure is used to define a Boundary condition. BCs affect
    the velocity of particles in a body.
    Each BC consists of two parts: a condition and an effect.


    Condition: The condition is used to determine which particles (in a body)
    will have the BC applied to them. The condition consists of....
        Position_Type: Either Reference or Spatial. Defines if the BC is applied
        relative to the particle's reference or spatial positions.

        Normal_Plane_Vector: Normal vector of a plane in 3d space

        Distance_To_Plane: The distance, in the direction of the Normal_Vector,
        from the origin to the plane. The plane consists of all points in 3d
        space such that
            <x,y,z> dot Condition_Normal_Vector = Condition_Distance_To_Plane

        Inequality: Defines how the BC is applied relative to the plane
        (LE = Less than or equal, GE = Greater than or Equal).

    To understand how this works, suppose that
        Condition_Position_Type = Reference
        Condition_Plane_Normal_Vector = {1,2,3}
        Condition_Plane_Distance = 15
        Condition_Inequality = GE
    Then any particle P whose reference position P_X satisifies
        Dot_Product(P_X, {1,2,3}) / Magnitude({1,2,3}) >= 15
    will have the BC effect applied to it.


    Effect: This determines the effect that is applied to all particles (in the
    body) that satisify the condition. The Effect consists of...
        Effect_x,y,z: Determines if there is an effect in the x, y, or z
        component. For example, if Effect_x = true, then the x component of the
        Effect_Vector will be applied to the particles that satisify the
        condition. If Effect_y is false then the particles that meet the
        condition will NOT have a BC applied to their y components.

        Effect_Vector: This determines how the BC affects the velocity
        of the particles that satisify the condition. The components of this
        vector whose corresponding Effect_x/y/z variable is true are applied.
        Those whose corresponding Effect_x/y/z variable is false are ignored. */

    // Condition
    Position_Type            Condition_Position_Type;
    Vector                   Condition_Plane_Normal_Vector;
    double                   Condition_Plane_Distance;
    Inequality               Condition_Inequality;

    // Effect
    bool                     Effect_x,
                             Effect_y,
                             Effect_z;
    Vector                   Effect_Vector;
  }; // struct Boundary_Condition {

  void Set_General_BCs(Body & Body_In,                     // The body we're applying the BC to
                       Array<Boundary_Condition> & BCs_In);        // The BCs being applied
  void Set_Box_BCs(Body & Box,                             // Reference to the box body
                   Box_Properties & Box_Parameters);       // Box's parameters
  void Set_Box_Particle_BCs(Particle & P_In,               // Particle that we're applying the BC to
                            Vector BC);                    // The BC that's being applied
} // namespace Simulation {

#endif
