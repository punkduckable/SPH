#if !defined(PARTICLE_ARRAY_SOURCE)
#define PARTICLE_ARRAY_SOURCE

#include "Particle_Array.h"

// Constructors, destructor
Particle_Array::Particle_Array(void) {
  Num_Particles = 0;
  Array = NULL;
} // Particle_Array::Particle_Array(void) {

Particle_Array::Particle_Array(const unsigned int Num_Particles_In) {
  if(Num_Particles == 0) {
    printf("An array of particles must have AT LEAST 1 particle\n");
    return;
  } // if(Num_Particles == 0) {

  Num_Particles = Num_Particles_In;
  Array = new Particle[Num_Particles];

  // Now assign each particle's Id
  for(unsigned int i = 0; i < Num_Particles; i++)
    Array[i].Set_ID(i);
} // Particle_Array::Particle_Array(const unsigned int Num_Particles_In) {

Particle_Array::~Particle_Array(void) {
  // Only attempt to delete the Array if it has been setup.
  if(Num_Particles != 0)
    delete [] Array;
} // Particle_Array::~Particle_Array(void) {



// Set methods
void Particle_Array::Set_Num_Particles(const unsigned int Num_Particles_In) {
  if(Num_Particles_In == 0) {
    printf("An array of particles must have AT LEAST 1 particle\n");
    return;
  } // if(Num_Particles_In == 0) {

  if(Num_Particles != 0) {
    printf("This particle array has already been setup!!!\n");
    printf("You can't change the number of particles in a Particle_Array!!!\n");
    return;
  }

  Num_Particles = Num_Particles_In;
  Array = new Particle[Num_Particles];

  // Now assign each particle's Id
  for(unsigned int i = 0; i < Num_Particles; i++)
    Array[i].Set_ID(i);
} // void Particle_Array::Set_Num_Particles(const unsigned int Num_Particles_In) {

void Particle_Array::Set_Inter_Particle_Spacing(const double IPS) {
  if(IPS <= 0) {
    printf("The inter particle spacing must be positive!\n");
    return;
  } // if(IPS <= 0) {

  // If Support_Radius has been set (meaning it's non-zero) then we can also
  // set h.
  if(Support_Radius != 0)
    Set_h(IPS*Support_Radius);

  Inter_Particle_Spacing = IPS;
} // void Particle_Array::Set_Inter_Particle_Spacing(const double IPS) {

void Particle_Array::Set_Support_Radius(const unsigned int SR_In) {
  if(SR_In == 0){
    printf("The support radius must be non-zero!\n");
    return;
  }

  // If the Inter particle spacing has been set (meaning it's non-zero) then
  // we can also set h.
  if(Inter_Particle_Spacing != 0)
    Set_h(SR_In*Inter_Particle_Spacing);

  Support_Radius = SR_In;
} // void Particle_Array::Set_Support_Radius(const unsigned int SR_In) {

#endif
