#if !defined(TESTS_HEADER)
#define TESTS_HEADER

#define X_SIDE_LENGTH 10
#define Y_SIDE_LENGTH 21
#define Z_SIDE_LENGTH 10
#define SUPPORT_RADIUS 3
#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif

// If true, this directs the test code to read from the 'Particles_File' file rather than generate an array of particles
const unsigned char Read_From_Particles_File = 1;
const unsigned char Save_Data_To_File = 0;

////////////////////////////////////////////////////////////////////////////////
// Initialize static particle class members. Note, if Read_From_Particles_File
// is true (1) then these values will be overwritten. 

double Particle::Inter_Particle_Spacing = 1;                                   //        : mm
double Particle::h = SUPPORT_RADIUS*Inter_Particle_Spacing;    // Support function radius   : mm
double Particle::Shape_Function_Amp = 15./(PI*h*h*h*h*h*h);    // Shape function Amplitude  : 1/mm^3

double Particle::Lame = 1.125;                       // Lame parameter             : Mpa
double Particle::mu0 = .275;                         // Shear modulus              : Mpa

double Particle::mu = 5e-5;                          // Viscosity                  : Mpa*s

double Particle::E = 0.770982;                       // Youngs modulus/Hourglass stiffness   : Mpa
double Particle::alpha = 7.5;                        // Hg control parameter       : Unitless

double Particle::Tau = .15;                          // Damage rate parameter      : unitless

void Vector_Tests(void);
void Tensor_Tests(void);
void List_Tests(void);
void Particle_Tests(void);
void Timing_Tests(void);

#endif
