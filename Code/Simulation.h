#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif
#define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000.)                       // Conversion from cpu cycles to msec

////////////////////////////////////////////////////////////////////////////////
// Initialize static particle class members. Note, these are set in the
// simulation namespace.
double Particle::Inter_Particle_Spacing;                                       //        : mm

double Particle::h;                                        // Support Radius in mm       : mm
unsigned int Particle::Support_Radius;                     // Support Radius in Inter Particle Spacings      : unitless
double Particle::Shape_Function_Amp;                       // Shape function Amplitude   : 1/mm^3

double Particle::Lame;                                     // Lame parameter             : Mpa
double Particle::mu0;                                      // Shear modulus              : Mpa

double Particle::mu;                                       // Viscosity                  : Mpa*s

double Particle::E;                                        // Youngs modulus/Hourglass stiffness   : Mpa
double Particle::alpha;                                    // Hg control parameter       : Unitless

double Particle::Tau;                                      // Damage rate parameter      : unitless


// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  void SetUp_Body(Particle * Body,
                   const unsigned int Num_Particles_Body,
                   const double IPS);
  void SetUp_Boundary(Particle * Boundary,
                       const unsigned int Num_Particles_Boundary,
                       const double IPS);

  // Simulation flags/properties
  const unsigned char Load_Data_From_File = 0;
  const unsigned char Save_Data_To_File = 1;
  const unsigned char Print_Forces = 1;
  const unsigned char TimeSteps_Between_Prints = 100;

  // Simulation dimensions
  unsigned int X_SIDE_LENGTH = 10;
  unsigned int Y_SIDE_LENGTH = 10;
  unsigned int Z_SIDE_LENGTH = 10;

  // TimeStep paramters
  const double dt = .00001;                                // Time step                  : s
  const unsigned int Num_Steps = 10000;                    // Number of time steps


  // Default static particle class values (Note, these will be overwritten if
  // 'Read_From_Particles_File' is true (1))
  void Set_Static_Particle_Members(void) {
    const double IPS = 1;
    const unsigned int Support_Radius = 3;
    const double h = IPS*Support_Radius;

    Particle::Inter_Particle_Spacing = IPS;                                    //        : mm
    Particle::h = h;                                       // Support Radius in mm       : mm
    Particle::Support_Radius = Support_Radius;             // Support Radius in Inter Particle Spacings      : unitless
    Particle::Shape_Function_Amp = 15./(PI*h*h*h*h*h*h);   // Shape function Amplitude   : 1/mm^6

    Particle::Lame = 1.125;                                // Lame parameter             : Mpa
    Particle::mu0 = .275;                                  // Shear modulus              : Mpa

    Particle::mu = 5e-4;                                   // Viscosity                  : Mpa*s

    Particle::E = 0.770982;                                // Youngs modulus/Hourglass stiffness   : Mpa
    Particle::alpha = 7.5;                                 // Hg control parameter       : Unitless

    Particle::Tau = .15;                                   // Damage rate parameter      : unitless
  } //   void Set_Static_Particle_Members(void) {
} // namespace Simulation {



#endif
