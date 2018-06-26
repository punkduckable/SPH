#include <iostream>
#include <math.h>
#include <time.h>
#include <string>
#include <unistd.h>

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Definitions
#define PI 3.1415926535897932384626
#define T 2

// Header files
#include "Classes.h"
#include "Vector.h"
#include "Tensor.h"
#include "Particle.h"
#include "List.h"
#include "VTK_File.h"
#include "Tests.h"

// Prototypes
Tensor Dyadic_Product(const Vector & V1,const Vector & V2);

// Source files
#include "Vector.c"
#include "Tensor.c"
#include "Particle.c"
#include "List.c"
#include "VTK_File.c"
#include "Tests.c"

////////////////////////////////////////////////////////////////////////////////
// Initialize static members of particle class

double Particle::h = 2;                                    // Suppoer function raduius   : mm
double Particle::Shape_Function_Amp = 15./(PI*h*h*h*h*h*h);// Shape function amplitude   : mm^-3

double Particle::Lame = 9;                                 // Lame parameter             : Mpa
double Particle::mu0 = .1;                                 // Shear modulus              : Mpa

double Particle::mu = 5e-3;                                // Viscosity                  : Mpa*s

double Particle::E = .01;                                  // Hourglass stiffness        : Mpa
double Particle::alpha = 50;                               // Hg control parameter       : Unitless

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor S;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */

  /* Unrolled loop (runs slower than 1 rolled loop for some reason)
  T[3*0 + 0] = V1[0]*V2[0];            // i = 0, j = 0
  T[3*0 + 1] = V1[0]*V2[1];            // i = 0, j = 1
  T[3*0 + 2] = V1[0]*V2[2];            // i = 0, j = 2

  T[3*1 + 0] = V1[1]*V2[0];            // i = 1, j = 0
  T[3*1 + 1] = V1[1]*V2[1];            // i = 1, j = 1
  T[3*1 + 2] = V1[1]*V2[2];            // i = 1, j = 2

  T[3*2 + 0] = V1[2]*V2[0];            // i = 2, j = 0
  T[3*2 + 1] = V1[2]*V2[1];            // i = 2, j = 1
  T[3*2 + 2] = V1[2]*V2[2];            // i = 2, j = 2
  */

  // Old loop. (works better than 9 statements with O2 optimization)
  for(int i = 0; i < 3; i++) {
    S[3*i + 0] = V1[i]*V2[0];
    S[3*i + 1] = V1[i]*V2[1];
    S[3*i + 2] = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++)

  return S;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {
  // Run Vector, Tensor, List tests.
  //Vector_Tests();
  //Tensor_Tests();
  //List_Tests();
  Particle_Tests();
  //Timing_Tests();

  return 0;
} // int main() {
