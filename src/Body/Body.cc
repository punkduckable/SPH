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



void Body::Set_Boundary(const bool Boundary_In) { Is_Boundary = Boundary_In; }

void Body::Set_First_Time_Step(const bool First_In) { First_Time_Step = First_In; }

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

bool Body::Get_Cuboid(void) const { return Is_Cuboid; }
unsigned Body::Get_X_SIDE_LENGTH(void) const {
  assert( (*this).Is_Cuboid );
  return X_SIDE_LENGTH;
} // unsigned Get_X_SIDE_LENGTH(void) const {
unsigned Body::Get_Y_SIDE_LENGTH(void) const {
  assert( (*this).Is_Cuboid );
  return Y_SIDE_LENGTH;
} // unsigned Get_Y_SIDE_LENGTH(void) const {
unsigned Body::Get_Z_SIDE_LENGTH(void) const {
  assert( (*this).Is_Cuboid );
  return Z_SIDE_LENGTH;
} // unsigned Get_Z_SIDE_LENGTH(void) const {

bool Body::Get_Boundary(void) const { return Is_Boundary; }

bool Body::Get_First_Time_Step(void) const { return First_Time_Step; }





////////////////////////////////////////////////////////////////////////////////
// Printing Methods

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
  printf(         "Material:                     %s\n",    Body_Material.Name.c_str());
  printf(         "Lame:                         %lf\n",   Body_Material.Lame);
  printf(         "mu0 (Shear modulus):          %lf\n",   Body_Material.mu0);
  printf(         "mu (Viscosity):               %lf\n",   mu);
  printf(         "E (Young's modulus):          %lf\n",   Body_Material.E);
  printf(         "Tau (Damage rate):            %lf\n\n", Tau);
} // void Body::Print_Parameters(void) const {



void Body::Print_Net_External_Force(const unsigned time_step) {
  /* This function is used to find and print the net external force on a body
  This function can NOT be called by multiple threads at once (this
  function is not thread safe). */

  // First, open the file.
  FILE * File;
  if(Times_Printed_Net_External_Force == 0) {
    File = fopen("../Files/Force_Files/Net_External_Force.txt","w");
  } // if(Times_Printed_Net_External_Force == 0) {
  else {
    File = fopen("../Files/Force_Files/Net_External_Force.txt","a");
  } // else {

  // Increment the number of times that we're printed net force data.
  Times_Printed_Net_External_Force++;

  // Now add up net external force on supplied particle array and print it out
  // Note that we must do this using a single thread
  Vector Net_Contact_Force = {0,0,0};

  for(unsigned i = 0; i < Num_Particles; i++) {
    Net_Contact_Force += Particles[i].Get_Force_Friction();
    Net_Contact_Force += Particles[i].Get_Force_Contact();
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  fprintf(File,"%6u:  <%10.4f, %10.4f, %10.4f>\n", time_step, Net_Contact_Force(0), Net_Contact_Force(1), Net_Contact_Force(2));

  // Now close the file.
  fclose(File);
} // void Body::Print_Net_External_Force(const unsigned time_step) {



void Body::Print_Particle_Forces(void) {
  // Create a file path for the new file (based on the Body's name
  // and time_step)
  char Buf[6];
  sprintf(Buf,"%05u",Times_Printed_Particle_Forces);
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name.c_str();
  File_Path +=  "_Force_";
  File_Path +=  Buf;

  // Now open the file.
  FILE * File = fopen(File_Path.c_str(), "w");

  // Increment the number of times that we're printed particle force data.
  Times_Printed_Particle_Forces++;

  // Print header.
  fprintf(File,"  ID  |");
  fprintf(File," Particle Pos  |");
  fprintf(File,"        Internal Force        |");

  #if defined(PARTICLE_DEBUG)
    fprintf(File,"        Viscous Force         |");
  #endif

  fprintf(File,"        Contact Force         |");
  fprintf(File,"        Friction Force        |");
  fprintf(File,"        Hourglass Force       |");
  fprintf(File,"\n");

  // Cycle through particles, print spacial positions of each particle
  for(unsigned i = 0; i < Num_Particles; i++) {
    fprintf(File,"%6u|", Particles[i].Get_ID());
    fprintf(File,"%4.1f,%4.1f,%4.1f | ",    Particles[i].X[0],            Particles[i].X[1],            Particles[i].X[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Int[0],    Particles[i].Force_Int[1],    Particles[i].Force_Int[2]);

    #if defined(PARTICLE_DEBUG)
      fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Visc[0],   Particles[i].Force_Visc[1],   Particles[i].Force_Visc[2]);
    #endif

    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Contact[0], Particles[i].Force_Contact[1], Particles[i].Force_Contact[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Friction[0],Particles[i].Force_Friction[1],Particles[i].Force_Friction[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e>\n",   Particles[i].Force_HG[0],      Particles[i].Force_HG[1],      Particles[i].Force_HG[2]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Body::Print_Particle_Forces(void) {
