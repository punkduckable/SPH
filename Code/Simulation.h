#if !defined(SUMULATION_HEADER)
#define SUMULATION_HEADER

#define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000.)                       // Conversion from cpu cycles to msec

// Simulation napespace (stores the variables and functions that are needed to
// run a simulation)
namespace Simulation {
  void SetUp_Body(Particle_Array & Body);
  void SetUp_Boundary(Particle_Array & Boundary);
  // Simulation flags/properties
  const unsigned char Load_Data_From_File = 0;
  const unsigned char Save_Data_To_File = 1;
  const unsigned char Print_Forces = 1;
  const unsigned char TimeSteps_Between_Prints = 100;

  // Simulation dimensions
  unsigned int X_SIDE_LENGTH = 10;
  unsigned int Y_SIDE_LENGTH = 5;
  unsigned int Z_SIDE_LENGTH = 10;

  // Support radius
  unsigned int Support_Radius = 4;                         // Support radius in units of Inter Particle spacings

  // TimeStep paramters
  const double dt = .00001;                                // Time step                  : s
  const unsigned int Num_Steps = 10000;                    // Number of time steps

  // Default Particle_Array members
  void Set_Particle_Array_Members(Particle_Array & Particles) {
    const double IPS = 1;                                    // Default Inter particle spacing

    Particles.Set_Inter_Particle_Spacing(IPS);                                   //        : mm
    Particles.Set_Support_Radius(Support_Radius);            // Support Radius in Inter Particle Spacings      : unitless

    Particles.Set_Lame(1.125);                               // Lame parameter             : Mpa
    Particles.Set_mu0(.275);                                 // Shear modulus              : Mpa

    Particles.Set_mu(5e-4);                                  // Viscosity                  : Mpa*s

    Particles.Set_E(0.770982);                               // Youngs modulus/Hourglass stiffness   : Mpa
    Particles.Set_alpha(7.5);                                // Hg control parameter       : Unitless

    Particles.Set_Tau(.15);                                  // Damage rate parameter      : unitless
  } // void Set_Particle_Array_Members(Particle_Array & Particles) {
} // namespace Simulation {

#endif
