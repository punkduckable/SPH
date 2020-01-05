#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "IO/FEB_File.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

////////////////////////////////////////////////////////////////////////////////
// Define external simulation.h variables

// Body properties
namespace Simulation {
  unsigned Num_Bodies = 0;                       // Number of bodies in simulation
  std::string * Names = nullptr;                 // The names of each body (name must match File name if reading from FEB file)
  bool * Is_Cuboid = nullptr;                    // Which bodies are cuboids
  bool * Is_Boundary = nullptr;                  // Which bodies are boundaries (can be from FEB file or cuboid)
  bool * Is_Damagable = nullptr;                 // Which bodies can be damaged
  bool * From_FEB_File = nullptr;                // Which bodies will be read from file
  unsigned * Steps_Per_Update = nullptr;         // How many time steps pass between updating this Body's P-K tensor
  double * IPS = nullptr;                        // Inter particle spacing in mm.
  Vector * Dimensions = nullptr;                 // Dimensions of cuboids (only applicable for cuboids)
  Vector * Offset = nullptr;                     // Poisition offset (only applicable for cuboids)
  Vector * Initial_Velocity = nullptr;           // Initial velocity condition
  double * Num_Particles = nullptr;              // The number of particles in each body
  Materials::Material * Simulation_Materials = nullptr;    // Each bodies material
} // namespace Simulation {





////////////////////////////////////////////////////////////////////////////////
// Simulation set up

void Simulation::Body_Needle_Set_Up(void) {
  Num_Bodies                                   = 2;

  Names = new std::string[Num_Bodies];
  Is_Cuboid = new bool[Num_Bodies];
  Is_Boundary = new bool[Num_Bodies];
  Is_Damagable = new bool[Num_Bodies];
  From_FEB_File = new bool[Num_Bodies];
  Steps_Per_Update = new unsigned[Num_Bodies];
  IPS = new double[Num_Bodies];
  Dimensions = new Vector[Num_Bodies];
  Offset = new Vector[Num_Bodies];
  Initial_Velocity = new Vector[Num_Bodies];
  Num_Particles = new double[Num_Bodies];
  Simulation_Materials = new Materials::Material[Num_Bodies];

  Names[0]                                     = "Body";
  Is_Cuboid[0]                                 = true;
  Is_Boundary[0]                               = false;
  Is_Damagable[0]                              = true;
  From_FEB_File[0]                             = false;
  Steps_Per_Update[0]                          = 100;
  IPS[0]                                       = 1;
  Dimensions[0]                                = {20, 10, 20};
  Offset[0]                                    = {0, 0, 0};
  Initial_Velocity[0]                          = {0, 0, 0};
  Simulation_Materials[0]                      = Materials::Default;

  Names[1]                                     = "Needle";
  Is_Cuboid[1]                                 = true;
  Is_Boundary[1]                               = false;
  Is_Damagable[1]                              = false;
  From_FEB_File[1]                             = false;
  Steps_Per_Update[1]                          = 1;
  IPS[1]                                       = 1;
  Dimensions[1]                                = {4, 10, 4};
  Offset[1]                                    = {10-2, 11, 10-2};
  Initial_Velocity[1]                          = {0, -50, 0};
  Simulation_Materials[1]                      = Materials::Stainless_Steel;
} // void Simulation::Body_Needle_Set_Up(void) {



void Simulation::Set_Body_Members(Body & Body_In) {
  unsigned Support_Radius = 3;                   // Support radius in units of Inter Particle spacings

  Body_In.Set_Support_Radius(Support_Radius);    // Support Radius in Inter Particle Spacings      : unitless
  Body_In.Set_mu(1e-4);                          // Viscosity                  : Mpa*s
  Body_In.Set_alpha(.75);                        // Hg control parameter       : Unitless
  Body_In.Set_Tau(.15);                          // Damage rate parameter      : unitless
} // void Simulation::Set_Body_Members(Body & Body_In) {



void Simulation::Setup_Cuboid(Body & Body_In, const unsigned m) {
  unsigned i,j,k;
  TIME_TYPE time1;

  // Particle paramaters
  const double IPS = Body_In.Get_Inter_Particle_Spacing();                     //        : mm
  const double Particle_Volume = IPS*IPS*IPS;                                  //        : mm^3
  const double Particle_Radius = IPS*.578;                                     //        : mm
  const double Particle_Mass = Particle_Volume*Body_In.Get_density();          //        : g

  // Furst, let's get number of partilces in the Bodys
  const unsigned X_SIDE_LENGTH = Body_In.Get_X_SIDE_LENGTH();
  const unsigned Y_SIDE_LENGTH = Body_In.Get_Y_SIDE_LENGTH();
  const unsigned Z_SIDE_LENGTH = Body_In.Get_Z_SIDE_LENGTH();

  // Vectors to hold onto Parameters
  Vector X, x;
  Vector V = Initial_Velocity[m];                          // Initial_Velocity set in Simulation.h           : mm/s

  //////////////////////////////////////////////////////////////////////////////
  // Set up particles
  printf(         "\nGenerating particles for %s...",Body_In.Get_Name().c_str());
  time1 = Get_Time();

  // Set up Body_In
  /* Store particles in 'Vertical Column' major 'Row' semi-major order
  A vertical column is a set of particles with the same x and z coordinates,
  while a row is a set of particles with the same y and x coordinates. This
  ordering method places particles with the same */
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
        unsigned index = i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j;

        X = {i*IPS, j*IPS, k*IPS};
        X += Offset[m];
        x = X;                                                                 //        : mm

        Body_In[index].Set_Mass(Particle_Mass);                                //        : g
        Body_In[index].Set_Vol(Particle_Volume);                               //        : mm^3
        Body_In[index].Set_Radius(Particle_Radius);                            //        : mm
        Body_In[index].Set_X(X);                                               //        : mm
        Body_In[index].Set_x(x);                                               //        : mm
        Body_In[index].Set_V(V);                                               //        : mm/s
      } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  #if defined(_OPENMP)
    time1 = omp_get_wtime() - time1;
    printf(        "Done!\ntook %lf s\n",time1);
  #else
    time1 = clock() - time1;
    unsigned long MS_Gen = (unsigned long)(((float)time1)/((float)CLOCKS_PER_MS));
    printf(        "Done!\ntook %lu ms\n",MS_Gen);
  #endif

  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors (if the body is not a boundary)

  if(Body_In.Get_Boundary() == false) {
    printf(         "Generating %s's neighbor lists...", Body_In.Get_Name().c_str());
    time1 = Get_Time();
    Body_In.Find_Neighbors_Cuboid();

    #if defined(_OPENMP)
      time1 = omp_get_wtime() - time1;
      printf(        "Done!\ntook %lf s\n",time1);
    #else
      time1 = clock() - time1;
      unsigned long MS_Neighbor = (unsigned long)(((float)time1)/((float)CLOCKS_PER_MS));
      printf(       "Done!\ntook %lums\n",MS_Neighbor);
    #endif
  } //   if(Body_In.Get_Boundary() == false) {

  /*
  // Damage the 'cut'
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Body_In[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Body_In.Remove_Damaged_Particle(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2));
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {
  */
} // void Simulation::Setup_Cuboid(Body & Body_In, const unsigned m) {



void Simulation::Setup_FEB_Body(Body & FEB_Body, const unsigned m) {
  // First, we need to know how many particles we have, and the reference
  // position of each of the particles.
  Vector * X = nullptr;
  unsigned Num_Particles;
  FEB_File::Read_FEB_File(Names[m], &X, Num_Particles);    // Names in Simulation.h

  printf("\nReading in Particles for %s from FEB file...\n", FEB_Body.Get_Name().c_str());

  // Now we can set up the body
  FEB_Body.Set_Num_Particles(Num_Particles);

  // Now we can cycle through the particles, setting up each particle.
  const double IPS = FEB_Body.Get_Inter_Particle_Spacing();                    //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Mass = Particle_Volume*FEB_Body.Get_density();               //        : g

  Vector V = Initial_Velocity[m];                          // Initial_Velocity set in Simulation.h


  for(unsigned i = 0; i < Num_Particles; i++) {
    FEB_Body[i].Set_Mass(Particle_Mass);
    FEB_Body[i].Set_Vol(Particle_Volume);
    FEB_Body[i].Set_Radius(Particle_Radius);
    FEB_Body[i].Set_X(X[i]);
    FEB_Body[i].Set_x(X[i]);
    FEB_Body[i].Set_V(V);
  } //   for(unsigned i = 0; i < Num_Particles; i++) {

  // Now set up neighbors. (if the body is not a boundary)
  if(FEB_Body.Get_Boundary() == false) {
    printf("Setting up neighbors for %s...\n",FEB_Body.Get_Name().c_str());
    FEB_Body.Find_Neighbors();
    printf("Done!\n");
  } // if(FEB_Body.Get_Boundary() == false) {
} // void Simulation::Setup_FEB_Body(Body & FEB_Body, const unsigned m) {