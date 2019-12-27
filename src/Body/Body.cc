#include "Body.h"

// Set K (static member of Body class)
double Body::K = 400;

////////////////////////////////////////////////////////////////////////////////
// Constructors, destructor

Body::Body(void) {
  Array = NULL;
  Num_Particles = 0;
  X_SIDE_LENGTH = 0;
  Y_SIDE_LENGTH = 0;
  Z_SIDE_LENGTH = 0;
  Is_Cuboid = false;
  Is_Boundary = false;
  Inter_Particle_Spacing = 0;
  Support_Radius = 0;
  h = 0;
  Shape_Function_Amplitude = 0;
} // Body::Body(void) {



Body::Body(const unsigned Num_Particles_In) {
  if(Num_Particles == 0) {
    printf("An array of particles must have AT LEAST 1 particle\n");
    return;
  } // if(Num_Particles == 0) {

  Num_Particles = Num_Particles_In;
  Particles = new Particle[Num_Particles];

  // Now assign each particle's ID
  for(unsigned i = 0; i < Num_Particles; i++) { Particles[i].Set_ID(i); }

  // Set other members
  X_SIDE_LENGTH = 0;
  Y_SIDE_LENGTH = 0;
  Z_SIDE_LENGTH = 0;
  Is_Cuboid = false;
  Inter_Particle_Spacing = 0;
  Support_Radius = 0;
  h = 0;
  Shape_Function_Amplitude = 0;

  (*this).Particles_Set_Up = true;
} // Body::Body(const unsigned Num_Particles_In) {





////////////////////////////////////////////////////////////////////////////////
// Operator Overoading

Particle & Body::operator[](const unsigned i) {
  // Check that i is within bounds of Particles array.
  assert(i < (*this).Num_Particles);

  return (*this).Particles[i];
} // Particle & Body::operator[](const unsigned i) {



const Particle & Body::operator[](const unsigned i) const {
  // Check that i is within bounds of Particles Array.
  assert(i < (*this).Num_Particles);

  return (*this).Particles[i];
} // const Particle & Body::operator[](const unsigned i) const {





////////////////////////////////////////////////////////////////////////////////
// Set methods

void Body::Set_Num_Particles(const unsigned Num_Particles_In) {
  if(Is_Cuboid == true) {
    printf("This is a cuboid... set the x, y, and z side lengths\n");
    return;
  } // if(Is_Cuboid == true) {

  if(Num_Particles_In == 0) {
    printf("An array of particles must have AT LEAST 1 particle\n");
    return;
  } // if(Num_Particles_In == 0) {

  if((*this).Particles_Set_Up == true) {
    printf("This particle array has already been setup!!!\n");
    printf("You can't change the number of particles in a Body!!!\n");
    return;
  } // if((*this).Particles_Set_Up == true) {

  Num_Particles = Num_Particles_In;
  Particles = new Particle[Num_Particles];

  // Now assign each particle's Id
  for(unsigned i = 0; i < Num_Particles; i++) { Particles[i].Set_ID(i); }

  (*this).Particles_Set_Up = true;
} // void Body::Set_Num_Particles(const unsigned Num_Particles_In) {



void Body::Set_Inter_Particle_Spacing(const double IPS) {
  if(IPS <= 0) {
    printf("The inter particle spacing must be positive!\n");
    return;
  } // if(IPS <= 0) {

  // If Support_Radius has been set (meaning it's non-zero) then we can also
  // set h.
  if(Support_Radius != 0) { Set_h(IPS*Support_Radius); }

  Inter_Particle_Spacing = IPS;
} // void Body::Set_Inter_Particle_Spacing(const double IPS) {



void Body::Set_Support_Radius(const unsigned SR_In) {
  if(SR_In == 0) {
    printf("The support radius must be non-zero!\n");
    return;
  } // if(SR_In == 0){

  // If the Inter particle spacing has been set (meaning it's non-zero) then
  // we can also set h.
  if(Inter_Particle_Spacing != 0) { Set_h(SR_In*Inter_Particle_Spacing); }

  Support_Radius = SR_In;
} // void Body::Set_Support_Radius(const unsigned SR_In) {



void Body::Set_Cuboid_Dimensions(const Vector & Dimensions) {
  // Check if cuboid has already been set up
  if((*this).Particles_Set_Up == true) {
    printf("%s has already been set up! You can't change its dimensions\n", Name.c_str());
    return;
  } // if(Num_Particles != 0) {

  // Check for non-sensical input
  if(Dimensions(0) == 0 || Dimensions(1) == 0 || Dimensions(2) == 0) {
      printf("%s can't have a side length of zero...\n", Name.c_str());
      return;
  } // if(Dimensions(0) == 0 || Dimensions(1) == 0 || Dimensions(2) == 0) {

  // Now designate this Body as a Cuboid
  Is_Cuboid = true;

  X_SIDE_LENGTH = Dimensions(0);
  Y_SIDE_LENGTH = Dimensions(1);
  Z_SIDE_LENGTH = Dimensions(2);

  // Set up particles array.
  Num_Particles = X_SIDE_LENGTH*Y_SIDE_LENGTH*Z_SIDE_LENGTH;
  Particles = new Particle[Num_Particles];

  // Set ID, ijk coordinates of each particle in the newly allocated array.
  unsigned index = 0;
  for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {
    for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
      for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
        index++;
        Particles[index].Set_ID(index);
      } // for(unsigned k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(unsigned j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(unsigned i = 0; i < X_SIDE_LENGTH; i++) {

  (*this).Particles_Set_Up = true;
} // void Body::Set_Cuboid_Dimensions(const Vector & Dimensions); {



////////////////////////////////////////////////////////////////////////////////
// Other methods

void Body::Print_Parameters(void) const {
  printf(         "Name:                         %s\n",    Name.c_str());
  printf(         "Is a cuboid:                  %u\n",    (unsigned)Is_Cuboid);
  if(Is_Cuboid == true) {
    printf(       "X side length:                %u\n",    X_SIDE_LENGTH);
    printf(       "Y side length:                %u\n",    Y_SIDE_LENGTH);
    printf(       "Z side length:                %u\n",    Z_SIDE_LENGTH);
  } // if(Is_Cuboid) {

  printf(         "Number of particles:          %u\n",    Num_Particles);
  printf(         "Partciles Array address:      %p\n",    Particles);
  printf(         "Inter particle spacing:       %lf\n",   Inter_Particle_Spacing);
  printf(         "h:                            %lf\n",   h);
  printf(         "Support Radius:               %u\n",    Support_Radius);
  printf(         "Shape Function Amplitude:     %lf\n",   Shape_Function_Amplitude);
  printf(         "Material:                     %s\n",    Body_Matterial.Name.c_str());
  printf(         "Lame:                         %lf\n",   Body_Matterial.Lame);
  printf(         "mu0 (Shear modulus):          %lf\n",   Body_Matterial.mu0);
  printf(         "mu (Viscosity):               %lf\n",   mu);
  printf(         "E (Young's modulus):          %lf\n",   Body_Matterial.E);
  printf(         "Tau (Damage rate):            %lf\n\n", Tau);
} // void Body::Print_Parameters(void) const {
