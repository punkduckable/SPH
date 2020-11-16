#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "IO/FEB_File.h"
#include "Array.h"
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




// Declare Simulation Parameters
namespace Simulation {
  // IO paramaters
  bool Load_Simulation_From_Save;
  bool Save_Simulation_To_File;
  bool Print_Particle_Forces;
  bool Print_Body_Forces;
  unsigned TimeSteps_Between_Prints;

  // TimeStep paramters
  double dt;                                     // Time step        : s
  unsigned Num_Time_Steps;                       // Number of time steps

  // Contact
  double Contact_Distance;                       // Distance at which bodies begin contacting one another.   : mm
  double Friction_Coefficient;                                                 //        : unitless

  // Number of bodies
  unsigned Num_Bodies;                          // Number of bodies in simulation

  // Simulation setup parameters
  bool * From_FEB_File = nullptr;                // Which bodies will be read from file
  Array<General_Boundary_Condition> * General_BCs = nullptr;    // Specifies the general BCs for each body
  Vector * Position_Offset = nullptr;            // Position offset for particles in body
  Vector * Initial_Velocity = nullptr;           // Initial velocity condition
} // namespace Simulation {



void Simulation::Setup(Body ** Bodies) {
  /* Function description:
  This function, as the name implies, sets up a Simulation. This means setting
  all simulation parameters and creating the Bodies.
  Most of this is actually handeled by Load_Setup_File.

  Simulation::Run is the only thing that should call this function */

  TIME_TYPE time_load = Get_Time();
  printf("Setting up simulation...\n");

  // First, read in from the Setup file.
  *Bodies = Simulation::Load_Setup_File();

  /* Load Setup_File set "Load_Simulation_From_Save". If this parameter is true,
  then we should run Load_Simulation. */
  if(Simulation::Load_Simulation_From_Save == true) {
    IO::Load_Simulation(Bodies, Simulation::Num_Bodies);
  } // if(Simulation::Load_Simulation_From_Save == true) {

  /* Otherwise, we should finish setting up the bodies. */
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
                        i,(*Bodies)[i].Get_Name().c_str());
        throw Bad_Body_Setup(Buffer);
      } // if((*Bodies)[i].Get_Is_Box() == true && Simulation::From_FEB_File[i] == true) {

      // A body must either be a Box or be from file. If it's neither, then
      // we have no way of setting it up.
      if((*Bodies)[i].Get_Is_Box() == false && Simulation::From_FEB_File[i] == false) {
        char Buffer[500];
        sprintf(Buffer, "Bad Body Setup Exception: Thrown by Startup_Simulation\n"
                        "Body %d (named %s) is designated neither a Box nor from FEB file\n"
                        "However, each body must either be from FEB file or a Box (but not both)\n",
                        i,(*Bodies)[i].Get_Name().c_str());
        throw Bad_Body_Setup(Buffer);
      } //   if((*Bodies)[i].Get_Is_Box() == false && Simulation::From_FEB_File[i] == false) {



      //////////////////////////////////////////////////////////////////////////
      // Set up (*Bodies)[i]'s particles

      /* If the ith Body is a Box then set it up as a Box. Otherwise, if the ith
      both is from a FEB file, then read it in. */
      if((*Bodies)[i].Get_Is_Box() == true) { Simulation::Setup_Box((*Bodies)[i], i); }
      else {                                  Simulation::Setup_FEB_Body((*Bodies)[i], i); }



      //////////////////////////////////////////////////////////////////////////
      // Set up (*Bodies)[i]'s General Boundary Condition

      Simulation::Set_General_BCs((*Bodies)[i], General_BCs[i]);
    } // for(unsigned i = 0; i < Simulation::Num_Bodies; i++) {
  } // else {

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
  printf(         "Print_Body_Forces =           %u\n",    Simulation::Print_Body_Forces);
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
  delete [] General_BCs;
  delete [] Position_Offset;
  delete [] Initial_Velocity;
} // void Simulation::Setup(Body ** Bodies) {





void Simulation::Setup_Box(Body & Body_In, const unsigned m) {
  /* Function Description:

  This function sets up Body_In as a box using the parameters corresponding
  to body m in Bodies_Setup. This function should NOT be used if you're loading
  from a save.

  Simulation::Setup is the only thing that should call this function. */

  unsigned i,j,k;
  TIME_TYPE time1;

  // Particle paramaters
  const double IPS = Body_In.Get_Inter_Particle_Spacing();                     //        : mm
  const double Particle_Volume = IPS*IPS*IPS;                                  //        : mm^3
  const double Particle_Radius = IPS*.578;                                     //        : mm
  const double Particle_Mass = Particle_Volume*Body_In.Get_density();          //        : g

  // First, let's get number of partilces in Body_In
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

  Simulation::Setup is the only thing that should call this function */

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
