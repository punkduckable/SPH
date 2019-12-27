#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#define FRICTION_COEFFICIENT .1

#include "Body/Body.h"
#include "Materials.h"
#include "Vector/Vector.h"
#include "Tensor/Tensor.h"


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
  const unsigned TimeSteps_Between_Prints    = 1000;

  // TimeStep paramters
  const double dt                                = 5e-9;          // Time step        : s
  const unsigned Num_Steps                   = 100000;          // Number of time steps

  // Body properties
  unsigned Num_Arrays;                           // Number of bodies in simulation
  std::string * Names;                           // The names of each body (name must match File name if reading from FEB file)
  bool * Is_Cuboid;                              // Which bodies are cuboids
  bool * Is_Boundary;                            // Which bodies are boundaries (can be from FEB file or cuboid)
  bool * Is_Damagable;                           // Which bodies can be damaged
  bool * From_FEB_File;                          // Which bodies will be read from file
  unsigned * Steps_Per_Update;               // How many time steps pass between updating this Body's P-K tensor
  double * IPS;                                  // Inter particle spacing in mm.
  Vector * Dimensions;                           // Dimensions of cuboids (only applicable for cuboids)
  Vector * Offset;                               // Poisition offset (only applicable for cuboids)
  Vector * Initial_Velocity;                     // Initial velocity condition
  double * Num_Particles;                        // The number of particles in each body
  Materials::Material * Materials;               // Each bodies material

  // Contact parameter
  const double Contact_Distance = 1;             // Distance at which bodies begin contacting one another.   : mm

  void Use_Arrays_From_Code(void) {
    Num_Arrays                                   = 2;

    Names = new std::string[Num_Arrays];
    Is_Cuboid = new bool[Num_Arrays];
    Is_Boundary = new bool[Num_Arrays];
    Is_Damagable = new bool[Num_Arrays];
    From_FEB_File = new bool[Num_Arrays];
    Steps_Per_Update = new unsigned[Num_Arrays];
    IPS = new double[Num_Arrays];
    Dimensions = new Vector[Num_Arrays];
    Offset = new Vector[Num_Arrays];
    Initial_Velocity = new Vector[Num_Arrays];
    Num_Particles = new double[Num_Arrays];
    Materials = new Materials::Material[Num_Arrays];

    Names[0]                                     = "Body";
    Is_Cuboid[0]                                 = true;
    Is_Boundary[0]                               = false;
    Is_Damagable[0]                              = true;
    From_FEB_File[0]                             = false;
    Steps_Per_Update[0]                          = 100;
    IPS[0]                                       = 1;
    Dimensions[0]                                = {20, 10, 20};
    Offset[0]                                    = {0, 0, 0};
    Initial_Velocity[0]                          = {0, 0, 0};
    Materials[0]                                 = Materials::Default;

    Names[1]                                     = "Needle";
    Is_Cuboid[1]                                 = true;
    Is_Boundary[1]                               = false;
    Is_Damagable[1]                              = false;
    From_FEB_File[1]                             = false;
    Steps_Per_Update[1]                          = 1;
    IPS[1]                                       = 1;
    Dimensions[1]                                = {4, 10, 4};
    Offset[1]                                    = {10-2, 11, 10-2};
    Initial_Velocity[1]                          = {0, -50, 0};
    Materials[1]                                 = Materials::Stainless_Steel;
  } // void Use_Arrays_From_Code(void) {


  // Default Body members
  void Set_Body_Members(Body & Particles);
} // namespace Simulation {

#endif
