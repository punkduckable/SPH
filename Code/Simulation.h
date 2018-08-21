#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#define FRICTION_COEFFICIENT .25

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  void Setup_Cuboid(Particle_Array & Body, const unsigned int m);
  void Setup_FEB_Body(Particle_Array & FEB_Body, const unsigned int m);

  // Simulation flags/properties
  const unsigned char Load_Data_From_File        = 0;
  const unsigned char Save_Data_To_File          = 1;
  const unsigned char Print_Forces               = 0;
  const unsigned char Print_Net_Force            = 1;
  const unsigned int TimeSteps_Between_Prints    = 1000;

  // TimeStep paramters
  const double dt                                = .000001;// Time step                  : s
  const unsigned int Num_Steps                   = 100000; // Number of time steps

  // Particle_Array properties
  unsigned Num_Arrays;                           // Number of bodies in simulation
  std::string * Names;                           // The names of each body (name must match File name if reading from FEB file)
  bool * Is_Cuboid;                              // Which bodies are cuboids
  bool * Is_Boundary;                            // Which bodies are boundaries (can be from FEB file or cuboid)
  bool * Is_Damagable;                           // Which bodies can be damaged
  bool * From_FEB_File;                          // Which bodies will be read from file
  unsigned int * Steps_Between_Update;           // How many time steps pass between updating this Body's P-K tensor
  Vector * Dimensions;                           // Dimensions of cuboids (only applicable for cuboids)
  Vector * Offset;                               // Poisition offset (only applicable for cuboids)
  Vector * Initial_Velocity;                     // Initial velocity condition
  double * Num_Particles;                        // The number of particles in each body
  Materials::Material * Materials;               // Each bodies material

  void Use_Arrays_From_Code(void) {
    Num_Arrays                                   = 2;

    Names = new std::string[Num_Arrays];
    Is_Cuboid = new bool[Num_Arrays];
    Is_Boundary = new bool[Num_Arrays];
    Is_Damagable = new bool[Num_Arrays];
    From_FEB_File = new bool[Num_Arrays];
    Steps_Between_Update = new unsigned int[Num_Arrays];
    Dimensions = new Vector[Num_Arrays];
    Offset = new Vector[Num_Arrays];
    Initial_Velocity = new Vector[Num_Arrays];
    Num_Particles = new double[Num_Arrays];
    Materials = new Materials::Material[Num_Arrays];

    Names[0]                                     = "Body";
    Is_Cuboid[0]                                 = true;
    Is_Boundary[0]                               = false;
    Is_Damagable[0]                              = false;
    From_FEB_File[0]                             = false;
    Steps_Between_Update[0]                      = 1;
    Dimensions[0]                                = {10, 10, 10};
    Offset[0]                                    = {0, 10, 0};
    Initial_Velocity[0]                          = {0, 0, 0};
    Materials[0]                                 = Materials::Stainless_Steel;

    Names[1]                                     = "Boundary";
    Is_Cuboid[1]                                 = true;
    Is_Boundary[1]                               = true;
    Is_Damagable[1]                              = false;
    From_FEB_File[1]                             = false;
    Steps_Between_Update[1]                      = 1;
    Dimensions[1]                                = {15, 1, 20};
    Offset[1]                                    = {7, 0, -5};
    Initial_Velocity[1]                          = {0, 0, 0};
    Materials[1]                                 = Materials::Default;
  } // void Use_Arrays_From_Code(void) {


  // Default Particle_Array members
  void Set_Particle_Array_Members(Particle_Array & Particles) {
    unsigned int Support_Radius = 3;                       // Support radius in units of Inter Particle spacings
    double IPS = 1;                                        // Inter particle spacing     : mm

    Particles.Set_Inter_Particle_Spacing(IPS);                                 //        : mm
    Particles.Set_Support_Radius(Support_Radius);          // Support Radius in Inter Particle Spacings      : unitless

    Particles.Set_mu(1e-4);                                // Viscosity                  : Mpa*s

    Particles.Set_alpha(.75);                               // Hg control parameter       : Unitless

    Particles.Set_Tau(.15);                                // Damage rate parameter      : unitless
  } // void Set_Particle_Array_Members(Particle_Array & Particles) {
} // namespace Simulation {

#endif
