#include <stdio.h>
#include <math.h>
#include <unistd.h>

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Definitions
#define PI 3.1415926535897932384626

// Header files
#include "Vector.h"
#include "Tensor.h"
#include "Particle.h"
#include "List.h"
#include "Tests.h"

// Prototypes
Tensor Dyadic_Product(const Vector & V1,const Vector & V2);

// Source files
#include "Vector.c"
#include "Tensor.c"
#include "Particle.c"
#include "List.c"
#include "Tests.c"

////////////////////////////////////////////////////////////////////////////////
// Initialize static members of particle class

double Particle::density = 1;
double Particle::E = 1;
double Particle::alpha = 50;                    // What was used in the paper
Vector M_Set;
Vector Particle::M = M_Set;
double Particle::mu = 1;
double Particle::mu0 = 1;
double Particle::k1 = 1;
double Particle::k2 = 1;
double h = 5;
double Particle::h = h;
double Particle::A = 15./(PI*h*h*h*h*h*h);

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor T;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */
  for(int i = 0; i < 3; i++) {
    T(i,0) = V1[i]*V2[0];
    T(i,1) = V1[i]*V2[1];
    T(i,2) = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++) {

  return T;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {
  // Run Vector, Tensor, List tests.
  //Vector_Tests();
  //Tensor_Tests();
  //List_Tests();

  // Declare an array of particles
  Particle Particles[27];
  Particle & Current_Particle = Particles[0];

  // Initialize particle masses, volumes, etc..
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        Vector X = {(double)i,(double)j,(double)k};
        Vector x = X;
        double Mass = 1;
        double Vol = 1;

        Vector vel;
        if(i == 1 && j == 1 && k == 0)
          vel = {0.,0.,5.};
        else
          vel = {0.,0.,0.};

        Current_Particle = Particles[9*i + 3*j + k];
        Current_Particle.Set_Mass(Mass);
        Current_Particle.Set_Vol(Vol);
        Current_Particle.Set_X(X);
        Current_Particle.Set_x(x);
        Current_Particle.Set_vel(vel);
        sleep(1);
      }
    }
  }

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        Print(Particles[9*i + 3*j + k]);
      }
    }
  }
  return 0;
} // int main() {
