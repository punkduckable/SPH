#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000.)                       // Conversion from cpu cycles to msec

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  void Setup_Body(Particle_Array & Body);
  void Setup_Boundary(Particle_Array & Boundary);
  void Setup_FEB_Body(Particle_Array & FEB_Body, const std::string & File_Name);

  // Simulation flags/properties
  const unsigned char Load_Data_From_File = 0;
  const unsigned char Save_Data_To_File = 1;
  const unsigned char Print_Forces = 0;
  const unsigned char TimeSteps_Between_Prints = 100;

  // TimeStep paramters
  const double dt = .00001;                                // Time step                  : s
  const unsigned int Num_Steps = 0;                    // Number of time steps

  // Particle_Array properties
  unsigned Num_Arrays;
  std::string * Names;
  bool * Is_Cuboid;
  bool * Is_Boundary;
  bool * From_FEB_File;
  Vector * Dimensions;

  void Use_Arrays_From_Code(void) {
    Num_Arrays                                   = 2;
    Names = new std::string[Num_Arrays];
    Is_Cuboid = new bool[Num_Arrays];
    Is_Boundary = new bool[Num_Arrays];
    From_FEB_File = new bool[Num_Arrays];
    Dimensions = new Vector[Num_Arrays];

    Names[0]                                     = "Body";
    Is_Cuboid[0]                                 = true;
    Is_Boundary[0]                               = false;
    From_FEB_File[0]                             = false;
    Dimensions[0]                                = {20, 10, 20};

    Names[1]                                     = "Needle";
    Is_Cuboid[1]                                 = false;
    Is_Boundary[1]                               = false;
    From_FEB_File[1]                             = true;
    Dimensions[1]                                = {0,0,0};
  } // void Use_Arrays_From_Code(void) {


  // Default Particle_Array members
  void Set_Particle_Array_Members(Particle_Array & Particles) {
    unsigned int Support_Radius = 3;                         // Support radius in units of Inter Particle spacings
    double IPS = 1;                                          // Inter particle spacing     : mm

    Particles.Set_Inter_Particle_Spacing(IPS);                                 //        : mm
    Particles.Set_Support_Radius(Support_Radius);          // Support Radius in Inter Particle Spacings      : unitless

    Particles.Set_Lame(1.125);                             // Lame parameter             : Mpa
    Particles.Set_mu0(.275);                               // Shear modulus              : Mpa

    Particles.Set_mu(5e-4);                                // Viscosity                  : Mpa*s

    Particles.Set_E(0.770982);                             // Youngs modulus/Hourglass stiffness   : Mpa
    Particles.Set_alpha(7.5);                              // Hg control parameter       : Unitless

    Particles.Set_Tau(.15);                                // Damage rate parameter      : unitless
  } // void Set_Particle_Array_Members(Particle_Array & Particles) {
} // namespace Simulation {

#endif
