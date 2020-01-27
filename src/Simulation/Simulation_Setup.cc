#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "IO/FEB_File.h"
#include "Errors.h"
#include <assert.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif

// Prototypes for functions that are local to this file
namespace Simulation {
  void Setup_Box(Body & Body_In, const unsigned m);
  void Setup_FEB_Body(Body & FEB_Body, const unsigned m);
  void Bodies_Setup(void);                                 // Set up Body/Needle simulation
  void Set_Body_Members(Body & Body_In);                   // Set default body members
} // namespace Simulation {



void Simulation::Setup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {
  /* Function description:
  This function, as the name implies, Sets up a Simulation. This means setting
  all simulation parameters, the Bodies, and the Time_Step_Index array.
  Most of this is actually handeled by Load_Setup_File.

  In general, Run_Simulation is the only thing that should call this function */

  TIME_TYPE time_load = Get_Time();
  printf("Setting up simulation...\n");


  // First, read in from the Setup file.
  *Bodies = Simulation::Load_Setup_File();

  /* Load Setup_File set "Load_Simulation_From_Save". If this parameter is true,
  then we should run Load_Simulation. Otherwise, we should finish setting up
  the bodies. */
  if(Load_Simulation_From_Save == true) { IO::Load_Simulation(Bodies, Simulation::Num_Bodies); }

  else {
    for(unsigned i = 0; i < Simulation::Num_Bodies; i++) {
      //////////////////////////////////////////////////////////////////////////
      // Check for bad inputs!

      // A body can't both be a Box and be from an FEB file.
      if((*Bodies)[i].Get_Is_Box() == true && Simulation::From_FEB_File[i] == true) {
        char Buffer[500];
        sprintf(Buffer, "Bad Body Setup Exception: Thrown by Startup_Simulation\n"
                        "Body %d (named %s) is designated as both a Box and from FEB file\n"
                        "However, each body must either be from a FEB file or a Box (not both)\n",
                        i,Names[i].c_str());
        throw Bad_Body_Setup(Buffer);
      } // if((*Bodies)[i].Get_Is_Box() == true && Simulation::From_FEB_File[i] == true) {

      // A body must either be a Box or be from file. If it's neither, then
      // we have no way of setting it up.
      if((*Bodies)[i].Get_Is_Box() == false && Simulation::From_FEB_File[i] == false) {
        char Buffer[500];
        sprintf(Buffer, "Bad Body Setup Exception: Thrown by Startup_Simulation\n"
                        "Body %d (named %s) is designated neither a Box nor from FEB file\n"
                        "However, each body must either be from FEB file or a Box (but not both)\n",
                        i,Names[i].c_str());
        throw Bad_Body_Setup(Buffer);
      } //   if((*Bodies)[i].Get_Is_Box() == false && Simulation::From_FEB_File[i] == false) {



      //////////////////////////////////////////////////////////////////////////
      // Set up (*Bodies)[i]'s particles

      /* If the ith Body is a Box then set it up as a Box. Otherwise, if the ith
      both is from a FEB file, then read it in. */
      if((*Bodies)[i].Get_Is_Box() == true) { Setup_Box((*Bodies)[i], i); }
      else { Setup_FEB_Body((*Bodies)[i], i); }
    } // for(unsigned i = 0; i < Simulation::Num_Bodies; i++) {
  } // else {


  /* Set up the Time_Step_Index array (Note: if Load_Simulation_From_Save = true
  then the number of bodies is not known until Load_Simulation is done) */
  *Time_Step_Index = new unsigned[Num_Bodies];
  for(unsigned i = 0; i < Num_Bodies; i++) {  (*Time_Step_Index)[i] = 0; }


  // Report setup time.
  time_load = Time_Since(time_load);
  #if defined(_OPENMP)
    printf(     "\nDone! Setting up the simulation took %lf s\n", time_load);
  #else
    unsigned long MS_Load = (unsigned long)((double)time_load / (double)CLOCKS_PER_MS);
    printf(       "\nDone! Setting up the simulation took %lu ms\n", MS_Load);
  #endif


  // Display that the simulation has begun
  printf(         "\nRunning a Simulation...\n");
  printf(         "Load_Simulation_From_Save =   %u\n",    Simulation::Load_Simulation_From_Save);
  printf(         "Save_Simulation_To_File =     %u\n",    Simulation::Save_Simulation_To_File);
  printf(         "Print_Particle_Forces =       %u\n",    Simulation::Print_Particle_Forces);
  printf(         "Print_Net_External_Forces =   %u\n",    Simulation::Print_Net_External_Forces);
  printf(         "TimeSteps_Between_Prints =    %u\n",    Simulation::TimeSteps_Between_Prints);
  printf(         "Parallel execution =          ");
  #if defined(_OPENMP)
    printf(       "true\n");
    printf(       "Number of procs =             %u\n",omp_get_num_procs());
  #else
    printf(       "false\n");
  #endif

  printf(         "\nRunning a simulation with the following %d bodies:\n", Simulation::Num_Bodies);
  for(unsigned i = 0; i < Num_Bodies; i++) { (*Bodies)[i].Print_Parameters(); }

  // Now that the simulation has been set up, we can free the Simulation
  // parameters
  delete [] From_FEB_File;
  delete [] Box_Boundary_Conditions;
  delete [] Position_Offset;
  delete [] Initial_Velocity;
} // void Simulation::Setup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {



////////////////////////////////////////////////////////////////////////////////
// Functions to set up the simulation

void Simulation::Bodies_Setup(void) {
  Num_Bodies                                   = 2;

  Names = new std::string[Num_Bodies];
  From_FEB_File = new bool[Num_Bodies];
  Is_Box = new bool[Num_Bodies];
  Box_Parameters = new Box_Properties[Num_Bodies];
  Is_Fixed = new bool[Num_Bodies];
  Is_Damagable = new bool[Num_Bodies];
  Time_Steps_Per_Update = new unsigned[Num_Bodies];
  IPS = new double[Num_Bodies];
  Position_Offset = new Vector[Num_Bodies];
  Initial_Velocity = new Vector[Num_Bodies];
  Simulation_Materials = new Materials::Material[Num_Bodies];

  Names[0]                                     = "Body";
  Is_Box[0]                                    = true;
  Is_Fixed[0]                                  = false;
  Is_Damagable[0]                              = true;
  From_FEB_File[0]                             = false;
  Time_Steps_Per_Update[0]                     = 10;
  IPS[0]                                       = 1;
  Box_Parameters[0].Dimensions                 = {10, 5, 10};
  Box_Parameters[0].x_plus_BC                  = {0, FREE, FREE};
  Box_Parameters[0].x_minus_BC                 = {0, FREE, FREE};
  Box_Parameters[0].y_plus_BC                  = {FREE, 0, FREE};
  Box_Parameters[0].y_minus_BC                 = {FREE, 0, FREE};
  Box_Parameters[0].z_plus_BC                  = {FREE, FREE, 0};
  Box_Parameters[0].z_minus_BC                 = {FREE, FREE, 0};
  Position_Offset[0]                           = {0,0,0};
  Initial_Velocity[0]                          = {0, 0, 0};
  Simulation_Materials[0]                      = Materials::Default;

  Names[1]                                     = "Needle";
  Is_Box[1]                                    = true;
  Is_Fixed[1]                                  = false;
  Is_Damagable[1]                              = false;
  From_FEB_File[1]                             = false;
  Time_Steps_Per_Update[1]                     = 1;
  IPS[1]                                       = 1;
  Box_Parameters[1].Dimensions                 = {4, 10, 4};
  Box_Parameters[1].x_plus_BC                  = {FREE, FREE, FREE};
  Box_Parameters[1].x_minus_BC                 = {FREE, FREE, FREE};
  Box_Parameters[1].y_plus_BC                  = {0, -50, 0};
  Box_Parameters[1].y_minus_BC                 = {FREE, FREE, FREE};
  Box_Parameters[1].z_plus_BC                  = {FREE, FREE, FREE};
  Box_Parameters[1].z_minus_BC                 = {FREE, FREE, FREE};
  Position_Offset[1]                           = {5-2, 5.1, 5-2};
  Initial_Velocity[1]                          = {0, -500, 0};
  Simulation_Materials[1]                      = Materials::Old_Needle;
} // void Simulation::Bodies_Setup(void) {



void Simulation::Set_Body_Members(Body & Body_In) {
  unsigned Support_Radius = 3;                   // Support radius in units of Inter Particle spacings

  Body_In.Set_Support_Radius(Support_Radius);    // Support Radius in Inter Particle Spacings      : unitless
  Body_In.Set_mu(1e-4);                          // Viscosity                  : Mpa*s
  Body_In.Set_alpha(.75);                        // Hg control parameter       : Unitless
  Body_In.Set_Tau(.15);                          // Damage rate parameter      : unitless
} // void Simulation::Set_Body_Members(Body & Body_In) {



void Simulation::Setup_Box(Body & Body_In, const unsigned m) {
  /* Function Description:

  This function sets up Body_In as a box using the parameters corresponding
  to body m in Bodies_Setup. This function should NOT be used if you're loading
  from a save.

  In general, Setup_Simulation is the only thing that should call this function. */

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
  Vector V = Simulation::Initial_Velocity[m];              // Initial_Velocity set in Simulation.h           : mm/s

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
        X += Simulation::Position_Offset[m];
        x = X;                                                                 //        : mm

        Body_In[index].Set_Mass(Particle_Mass);                                //        : g
        Body_In[index].Set_Volume(Particle_Volume);                            //        : mm^3
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
  // Set up BCs
  Set_Box_BCs(Body_In, Simulation::Box_Boundary_Conditions[m]);



  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors (if the body is not fixed in place)
  if(Body_In.Get_Is_Fixed() == false) {
    printf(         "Generating %s's neighbor lists...", Body_In.Get_Name().c_str());
    time1 = Get_Time();
    Body_In.Find_Neighbors_Box();

    time1 = Time_Since(time1);
    #if defined(_OPENMP)
      printf(        "Done!\ntook %lf s\n",time1);
    #else
      unsigned long MS_Neighbor = (unsigned long)(((float)time1)/((float)CLOCKS_PER_MS));
      printf(       "Done!\ntook %lums\n",MS_Neighbor);
    #endif
  } //   if(Body_In.Get_Is_Fixed() == false) {

  /*
  // Damage the 'cut'
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Body_In[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Body_In.Remove_Damaged_Particle(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2));
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {
  */
} // void Simulation::Setup_Box(Body & Body_In, const unsigned m) {



void Simulation::Setup_FEB_Body(Body & FEB_Body, const unsigned m) {
  /* Function description:
  This function sets up a Body by reading in information from a FEB file.

  In general, Setup_Simulation is the only thing that should call this function */

  printf("\nReading in Particles for %s from FEB file...\n", FEB_Body.Get_Name().c_str());
  assert(Simulation::From_FEB_File[m] == true);


  // First, we need to know how many particles we have, and the reference
  // position of each of the particles.
  Vector * X = nullptr;
  unsigned Num_Particles;
  IO::Read_FEB_File(FEB_Body.Get_Name(), &X, Num_Particles);    // Names in Simulation.h


  // Now we can set up the body
  FEB_Body.Set_Num_Particles(Num_Particles);


  // Now we can cycle through the particles, setting up each particle.
  const double IPS = FEB_Body.Get_Inter_Particle_Spacing();                    //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Mass = Particle_Volume*FEB_Body.Get_density();               //        : g
  Vector V = Simulation::Initial_Velocity[m];              // Initial_Velocity set in Simulation.h


  for(unsigned i = 0; i < Num_Particles; i++) {
    FEB_Body[i].Set_Mass(Particle_Mass);
    FEB_Body[i].Set_Volume(Particle_Volume);
    FEB_Body[i].Set_Radius(Particle_Radius);
    FEB_Body[i].Set_X(X[i] + Simulation::Position_Offset[m]);
    FEB_Body[i].Set_x(X[i] + Simulation::Position_Offset[m]);
    FEB_Body[i].Set_V(V);
  } //   for(unsigned i = 0; i < Num_Particles; i++) {


  // Now set up neighbors. (if the body is not fixed in place)
  if(FEB_Body.Get_Is_Fixed() == false) {
    printf("Setting up neighbors for %s...\n",FEB_Body.Get_Name().c_str());
    FEB_Body.Find_Neighbors();
    printf("Done!\n");
  } // if(FEB_Body.Get_Is_Fixed() == false) {

  // Now free X.
  delete [] X;
} // void Simulation::Setup_FEB_Body(Body & FEB_Body, const unsigned m) {
