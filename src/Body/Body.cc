#include "Body.h"

// Set K (static member of Body class)
double Body::K = 400;

////////////////////////////////////////////////////////////////////////////////
// Constructors, destructor

Body::Body(void) {
  Array = nullptr;
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
  printf(         "Material:                     %s\n",    Body_Matterial.Name.c_str());
  printf(         "Lame:                         %lf\n",   Body_Matterial.Lame);
  printf(         "mu0 (Shear modulus):          %lf\n",   Body_Matterial.mu0);
  printf(         "mu (Viscosity):               %lf\n",   mu);
  printf(         "E (Young's modulus):          %lf\n",   Body_Matterial.E);
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



void Body::Export_Particle_Forces(const unsigned time_step) const {
  // Create a file path for the new file (based on the Body's name
  // and time_step)
  char Buf[6];
  sprintf(Buf,"%06u",time_step);
  std::string File_Path = "./IO/Force_Files/";
  File_Path += (*this).Name.c_str();
  File_Path +=  "_Force_";
  File_Path +=  Buf;

  // Now open the file.
  FILE * File = fopen(File_Path.c_str(), "w");

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
} // void Body::Export_Particle_Forces(void) const {
