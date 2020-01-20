#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "IO/FEB_File.h"
#include "Errors.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

////////////////////////////////////////////////////////////////////////////////
// Define external simulation.h variables

// Body properties
namespace Simulation {
  unsigned Num_Bodies = 0;                       // Number of bodies in simulation
  std::string * Names = nullptr;                 // The names of each body (name must match File name if reading from FEB file)
  bool * Is_Box = nullptr;                       // Which bodies are Boxs
  Box_Properties * Box_Parameters = nullptr;     // Specifies the dimensions, and BCs of the box bodies.
  bool * Is_Boundary = nullptr;                  // Which bodies are boundaries (can be from FEB file or Box)
  bool * Is_Damagable = nullptr;                 // Which bodies can be damaged
  bool * From_FEB_File = nullptr;                // Which bodies will be read from file
  unsigned * Time_Steps_Between_Updates=nullptr; // How many time steps pass between updating this Body's P-K tensor
  double * IPS = nullptr;                        // Inter particle spacing in mm.
  Vector * Position_Offset = nullptr;            // Position offset for particles in body
  Vector * Initial_Velocity = nullptr;           // Initial velocity condition
  Materials::Material * Simulation_Materials = nullptr;    // Each bodies material
} // namespace Simulation {





////////////////////////////////////////////////////////////////////////////////
// Functions to set up the simulation

void Simulation::Bodies_Setup(void) {
  Num_Bodies                                   = 2;

  Names = new std::string[Num_Bodies];
  From_FEB_File = new bool[Num_Bodies];
  Is_Box = new bool[Num_Bodies];
  Box_Parameters = new Box_Properties[Num_Bodies];
  Is_Boundary = new bool[Num_Bodies];
  Is_Damagable = new bool[Num_Bodies];
  Time_Steps_Between_Updates = new unsigned[Num_Bodies];
  IPS = new double[Num_Bodies];
  Position_Offset = new Vector[Num_Bodies];
  Initial_Velocity = new Vector[Num_Bodies];
  Simulation_Materials = new Materials::Material[Num_Bodies];

  Names[0]                                     = "Body";
  Is_Box[0]                                    = true;
  Is_Boundary[0]                               = false;
  Is_Damagable[0]                              = true;
  From_FEB_File[0]                             = false;
  Time_Steps_Between_Updates[0]                = 10;
  IPS[0]                                       = 1;
  Box_Parameters[0].Dimensions                 = {20, 10, 20};
  Box_Parameters[0].x_plus_BC                  = {0, Free_BC_Box, Free_BC_Box};
  Box_Parameters[0].x_minus_BC                 = {0, Free_BC_Box, Free_BC_Box};
  Box_Parameters[0].y_plus_BC                  = {Free_BC_Box, 0, Free_BC_Box};
  Box_Parameters[0].y_minus_BC                 = {Free_BC_Box, 0, Free_BC_Box};
  Box_Parameters[0].z_plus_BC                  = {Free_BC_Box, Free_BC_Box, 0};
  Box_Parameters[0].z_minus_BC                 = {Free_BC_Box, Free_BC_Box, 0};
  Position_Offset[0]                           = {0,0,0};
  Initial_Velocity[0]                          = {0, 0, 0};
  Simulation_Materials[0]                      = Materials::Default;

  Names[1]                                     = "Needle";
  Is_Box[1]                                    = true;
  Is_Boundary[1]                               = false;
  Is_Damagable[1]                              = false;
  From_FEB_File[1]                             = false;
  Time_Steps_Between_Updates[1]                = 1;
  IPS[1]                                       = 1;
  Box_Parameters[1].Dimensions                 = {4, 10, 4};
  Box_Parameters[1].x_plus_BC                  = {Free_BC_Box, Free_BC_Box, Free_BC_Box};
  Box_Parameters[1].x_minus_BC                 = {Free_BC_Box, Free_BC_Box, Free_BC_Box};
  Box_Parameters[1].y_plus_BC                  = {0, -50, 0};
  Box_Parameters[1].y_minus_BC                 = {Free_BC_Box, Free_BC_Box, Free_BC_Box};
  Box_Parameters[1].z_plus_BC                  = {Free_BC_Box, Free_BC_Box, Free_BC_Box};
  Box_Parameters[1].z_minus_BC                 = {Free_BC_Box, Free_BC_Box, Free_BC_Box};
  Position_Offset[1]                           = {10-2, 10.01, 10-2};
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
        X += Position_Offset[m];
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
    Body_In.Find_Neighbors_Box();

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
} // void Simulation::Setup_Box(Body & Body_In, const unsigned m) {



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



void Simulation::Startup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {
  //Display that the simulation has begun
  printf(         "\nRunning a Simulation...\n");
  printf(         "Load_Simulation_From_Save =   %u\n",    Load_Simulation_From_Save);
  printf(         "Save_Simulation =             %u\n",    Save_Simulation);
  printf(         "TimeSteps_Between_Prints =    %u\n",    TimeSteps_Between_Prints);
  printf(         "Print_Particle_Forces =       %u\n",    Print_Particle_Forces);
  printf(         "Parallel execution =          ");
  #if defined(_OPENMP)
    printf(       "1\n");
    printf(       "Number of procs =             %u\n",omp_get_num_procs());
  #else
    printf(       "0\n");
  #endif

  // Are we running a new simulation or loading an existing one?
  if(Load_Simulation_From_Save == 1) {
    TIME_TYPE time1 = Get_Time();

    // If loading an existing simulation, read in bodies from file
    printf(       "\nLoading from file....");
    IO::Load_Simulation(Bodies, Num_Bodies);

    // Now set up the time step counters
    *Time_Step_Index = new unsigned[Num_Bodies];
    for(unsigned i = 0; i < Num_Bodies; i++) {
      (*Time_Step_Index)[i] = 0;
    } // for(unsigned i = 0; i < Num_Bodies; i++) {


    #if defined(_OPENMP)
      time1 = omp_get_wtime() - time1;
      printf(     "Done!\ntook %lf s\n", time1);
    #else
      time1 = clock() - time1;
      unsigned long MS_Load = (unsigned long)((double)time1 / (double)CLOCKS_PER_MS);
      printf(       "Done!\ntook %lu ms\n", MS_Load);
    #endif
  } //   if(Load_Simulation_From_Save == 1) {

  else if(Load_Simulation_From_Save == 0) {
    // Use Bodies defined in Simulation.h
    Bodies_Setup();

    // First, allocate the array of Bodys, time step counters
    *Bodies = new Body[Num_Bodies];
    *Time_Step_Index = new unsigned[Num_Bodies];
    for(unsigned i = 0; i < Num_Bodies; i++) {  (*Time_Step_Index)[i] = 0; }

    // Now set up each body using the paramaters in Simulation.h
    for(unsigned i = 0; i < Num_Bodies; i++) {
      // Set Partilce Body's name
      (*Bodies)[i].Set_Name(Names[i]);

      // Set inter particle spacing
      (*Bodies)[i].Set_Inter_Particle_Spacing(IPS[i]);

      // Now set the ith Body's material
      (*Bodies)[i].Set_Material(Simulation_Materials[i]);

      // Now set wheather or not the ith Body is damagable
      (*Bodies)[i].Set_Damageable(Is_Damagable[i]);

      // Set Time Steps Per Update
      (*Bodies)[i].Set_Time_Steps_Between_Updates(Time_Steps_Between_Updates[i]);

      // Now set other Body members.
      Set_Body_Members((*Bodies)[i]);



      //////////////////////////////////////////////////////////////////////////
      // Check for bad inputs!

      // A body can't both be a Box and be from an FEB file.
      if(Is_Box[i] == true && From_FEB_File[i] == true) {
        char Buffer[500];
        sprintf(Buffer, "Bad Body Setup Exception: Thrown by Startup_Simulation\n"
                        "Body %d (named %s) is designated as both a Box and from FEB file\n"
                        "However, each body must either be from a FEB file or a Box (not both)\n",
                        i,Names[i].c_str());
        throw Bad_Body_Setup(Buffer);
      } // if(Is_Box[i] == true && From_FEB_File[i] == true) {

      // A body must either be a Box or be from file. If it's neither, then
      // we have no way of setting it up.
      if(Is_Box[i] == false && From_FEB_File[i] == false) {
        char Buffer[500];
        sprintf(Buffer, "Bad Body Setup Exception: Thrown by Startup_Simulation\n"
                        "Body %d (named %s) is designated neither a Box nor from FEB file\n"
                        "However, each body must either be from FEB file or a Box (but not both)\n",
                        i,Names[i].c_str());
        throw Bad_Body_Setup(Buffer);
      } // if(Is_Box[i] == false && Is_Boundary == false) {



      //////////////////////////////////////////////////////////////////////////
      // Now set up the Body's particles

      /* If the body is a boundary, we need to designate it as such.
      note: this applies if the body is a Box, or from FEB file... anything
      can be a boundary! */
      if(Is_Boundary[i] == true) {
        (*Bodies)[i].Set_Boundary(true);
      } // if(Is_Boundary[i] == true) {

      // Set up body as a Box if it is a Box
      if(Is_Box[i] == true) {
        (*Bodies)[i].Set_Box_Dimensions(Box_Parameters[i].Dimensions[0],
                                        Box_Parameters[i].Dimensions[1],
                                        Box_Parameters[i].Dimensions[2]);
        Setup_Box((*Bodies)[i], i);
      } // if(Is_Box[i] == true) {

      // if the body is from file, read it in
      else if(From_FEB_File[i] == true) {
        Setup_FEB_Body((*Bodies)[i], i);
      } // else if(From_FEB_File[i] == true) {

    } // for(unsigned i = 0; i < Num_Bodies; i++) {
  } // else if(Load_Simulation_From_Save == 0) {

  // Now that the body is loaded, print paramaters.
  printf(         "\nRuinning with the following paramaters....\n");
  for(unsigned i = 0; i < Num_Bodies; i++) {
    (*Bodies)[i].Print_Parameters();
  } // for(unsigned i = 0; i < Num_Bodies; i++) {
} // void Simulation::Startup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {
