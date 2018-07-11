#if !defined(TESTS_HEADER)
#define TESTS_HEADER

#define X_SIDE_LENGTH 10
#define Y_SIDE_LENGTH 30
#define Z_SIDE_LENGTH 10
#define SUPPORT_RADIUS 3
#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif

////////////////////////////////////////////////////1////////////////////////////
// Initialize static particle class members

const double Particle::Inter_Particle_Spacing = 1;                            //        : mm
const double Particle::h = SUPPORT_RADIUS*Inter_Particle_Spacing;    // Support function radius   : mm
const double Particle::Shape_Function_Amp = 15./(PI*h*h*h*h*h*h);    // Shape function amplitude  : mm^-3

const double Particle::Lame = 1.125;                       // Lame parameter             : Mpa
const double Particle::mu0 = .275;                         // Shear modulus              : Mpa

const double Particle::mu = 5e-5;                          // Viscosity                  : Mpa*s

const double Particle::E = 0.770982;                       // Youngs modulus/Hourglass stiffness   : Mpa
const double Particle::alpha = 7.5;                        // Hg control parameter       : Unitless

const double Particle::Tau = .15;                          // Damage rate paramater      : unitless

void Vector_Tests(void);
void Tensor_Tests(void);
void List_Tests(void);
void Particle_Tests(void);
void Timing_Tests(void);

#endif
