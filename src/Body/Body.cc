#include "Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include <assert.h>
#include <math.h>

// Set K (static member of Body class)
double Body::K = 400;

////////////////////////////////////////////////////////////////////////////////
// Constructors, destructor

Body::Body(void) {
  Particles = nullptr;
  Num_Particles = 0;
  X_SIDE_LENGTH = 0;
  Y_SIDE_LENGTH = 0;
  Z_SIDE_LENGTH = 0;
  Is_Box = false;
  Is_Fixed = false;
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
  Is_Box = false;
  Inter_Particle_Spacing = 0;
  Support_Radius = 0;
  h = 0;
  Shape_Function_Amplitude = 0;

  (*this).Particles_Set_Up = true;
} // Body::Body(const unsigned Num_Particles_In) {



Body::~Body(void) { }





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
// Setters

//---//---//---//---//---//---//---//---//---//---//---//---//---//---//---//--//
// Private

void Body::Set_h(const double h_In) {
  (*this).h = h_In;
  (*this).Shape_Function_Amplitude =  15./(PI*pow(h_In,6));
} // void Body::Set_h(const double h_In) {

//---//---//---//---//---//---//---//---//---//---//---//---//---//---//---//--//
// Public

void Body::Set_Num_Particles(const unsigned Num_Particles_In) {
  if(Is_Box == true) {
    printf("This is a Box... set the x, y, and z side lengths\n");
    return;
  } // if(Is_Box == true) {

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



void Body::Set_Name(const std::string & S_In) { Name = S_In; }



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



void Body::Set_Material(const Materials::Material & Mat_In) { Body_Material = Mat_In; }
void Body::Set_mu(const double mu_In) { mu = mu_In; }
void Body::Set_alpha(const double alpha_In) { alpha = alpha_In; }

void Body::Set_Tau(const double Tau_In) { Tau = Tau_In; }
void Body::Set_Damageable(const bool D_In) { Damageable = D_In; }




void Body::Set_Box_Dimensions(const unsigned Dim_x, const unsigned Dim_y, const unsigned Dim_z) {
  // Check if Box has already been set up
  if((*this).Particles_Set_Up == true) {
    printf("%s has already been set up! You can't change its dimensions\n", Name.c_str());
    return;
  } // if(Num_Particles != 0) {

  // Check for non-sensical input
  if(Dim_x == 0 || Dim_y == 0 || Dim_z == 0) {
      printf("%s can't have a side length of zero...\n", Name.c_str());
      return;
  } // if(Dim_x == 0 || Dim_y == 0 || Dim_z == 0) {

  // Now designate this Body as a Box
  Is_Box = true;

  X_SIDE_LENGTH = Dim_x;
  Y_SIDE_LENGTH = Dim_y;
  Z_SIDE_LENGTH = Dim_z;

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
} // void Body::Set_Box_Dimensions(const unsigned Dim_x, const unsigned Dim_y, const unsigned Dim_z) {



void Body::Set_Is_Fixed(const bool Is_Fixed_In) { Is_Fixed = Is_Fixed_In; }

void Body::Set_First_Time_Step(const bool First_In) { First_Time_Step = First_In; }
void Body::Set_Time_Steps_Between_Updates(const unsigned Steps_In) { Time_Steps_Between_Updates = Steps_In; }

void Body::Set_F_Index(const unsigned char i) {
  assert(i <= 1);
  F_Index = i;
} // void Body::Set_F_Counter(const unsigned char i) {



void Body::Increment_F_Index(void) {
  if(F_Index == 0) { F_Index++; }
  else { F_Index = 0; }
} // void Body::Increment_F_Counter(void) {





////////////////////////////////////////////////////////////////////////////////
// Getters

unsigned Body::Get_Num_Particles(void) const { return Num_Particles; }
std::string Body::Get_Name(void) const { return Name; }

double Body::Get_Inter_Particle_Spacing(void) const { return Inter_Particle_Spacing; }
unsigned Body::Get_Support_Radius(void) const { return Support_Radius; }

double Body::Get_h(void) const { return h; }
double Body::Get_Shape_Function_Amplitude(void) const { return Shape_Function_Amplitude; }
Materials::Material Body::Get_Material(void) const { return Body_Material; }
double Body::Get_Lame(void) const { return Body_Material.Lame; }
double Body::Get_mu0(void) const { return Body_Material.mu0; }
double Body::Get_mu(void) const { return mu; }
double Body::Get_E(void) const { return Body_Material.E; }
double Body::Get_density(void) const { return Body_Material.density; }
double Body::Get_alpha(void) const { return alpha; }

unsigned char Body::Get_F_Index(void) const { return F_Index; }

double Body::Get_Tau(void) const { return Tau; }
bool Body::Get_Damagable(void) const { return Damageable; }

bool Body::Get_Is_Box(void) const { return Is_Box; }
unsigned Body::Get_X_SIDE_LENGTH(void) const {
  assert( (*this).Is_Box );
  return X_SIDE_LENGTH;
} // unsigned Get_X_SIDE_LENGTH(void) const {
unsigned Body::Get_Y_SIDE_LENGTH(void) const {
  assert( (*this).Is_Box );
  return Y_SIDE_LENGTH;
} // unsigned Get_Y_SIDE_LENGTH(void) const {
unsigned Body::Get_Z_SIDE_LENGTH(void) const {
  assert( (*this).Is_Box );
  return Z_SIDE_LENGTH;
} // unsigned Get_Z_SIDE_LENGTH(void) const {

bool Body::Get_Is_Fixed(void) const { return Is_Fixed; }

bool Body::Get_First_Time_Step(void) const { return First_Time_Step; }
unsigned Body::Get_Time_Steps_Between_Updates(void) const { return Time_Steps_Between_Updates; }
