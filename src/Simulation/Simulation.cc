#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "IO/VTK_File.h"
#include "IO/Data_Dump.h"
#include "IO/FEB_File.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

namespace Simulation {
    // Body properties
    unsigned Num_Arrays;                           // Number of bodies in simulation
    std::string * Names;                           // The names of each body (name must match File name if reading from FEB file)
    bool * Is_Cuboid;                              // Which bodies are cuboids
    bool * Is_Boundary;                            // Which bodies are boundaries (can be from FEB file or cuboid)
    bool * Is_Damagable;                           // Which bodies can be damaged
    bool * From_FEB_File;                          // Which bodies will be read from file
    unsigned * Steps_Per_Update;                   // How many time steps pass between updating this Body's P-K tensor
    double * IPS;                                  // Inter particle spacing in mm.
    Vector * Dimensions;                           // Dimensions of cuboids (only applicable for cuboids)
    Vector * Offset;                               // Poisition offset (only applicable for cuboids)
    Vector * Initial_Velocity;                     // Initial velocity condition
    double * Num_Particles;                        // The number of particles in each body
    Materials::Material * Materials;               // Each bodies material
} // namespace Simulation {



void Simulation::Use_Arrays_From_Code(void) {
  Num_Arrays                                   = 2;

  Names = new std::string[Num_Arrays];
  Is_Cuboid = new bool[Num_Arrays];
  Is_Boundary = new bool[Num_Arrays];
  Is_Damagable = new bool[Num_Arrays];
  From_FEB_File = new bool[Num_Arrays];
  Steps_Per_Update = new unsigned[Num_Arrays];
  IPS = new double[Num_Arrays];
  Dimensions = new Vector[Num_Arrays];
  Offset = new Vector[Num_Arrays];
  Initial_Velocity = new Vector[Num_Arrays];
  Num_Particles = new double[Num_Arrays];
  Materials = new Materials::Material[Num_Arrays];

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
  Materials[0]                                 = Materials::Default;

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
  Materials[1]                                 = Materials::Stainless_Steel;
} // void Simulation::Use_Arrays_From_Code(void) {





void Simulation::Set_Body_Members(Body & Body_In) {
  unsigned Support_Radius = 3;                           // Support radius in units of Inter Particle spacings

  Body_In.Set_Support_Radius(Support_Radius);          // Support Radius in Inter Particle Spacings      : unitless
  Body_In.Set_mu(1e-4);                                // Viscosity                  : Mpa*s
  Body_In.Set_alpha(.75);                              // Hg control parameter       : Unitless
  Body_In.Set_Tau(.15);                                // Damage rate parameter      : unitless
} // void Simulation::Set_Body_Members(Body & Body_In) {




void Simulation::Run_Simulation(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned i,j,k,l,m;

  // Computation time measurement variables
  clock_t timer1,
          timer2,
          update_BC_timer = 0,
          update_P_timer = 0,
          contact_timer = 0,
          update_x_timer = 0,
          Print_timer = 0;

  // Set up Bodys.
  Body * Arrays;                                           // Will point to the Particle Array's for this simulation
  unsigned * Time_Step_Index;                                    // Time step counters for each particle array

  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.

  //Display that the simulation has begun
  printf(         "\nRunning a simulation...\n");
  printf(         "Load_Data_From_File =         %u\n",    Load_Data_From_File);
  printf(         "Save_Data_To_File =           %u\n",    Save_Data_To_File);
  printf(         "TimeSteps_Between_Prints =    %u\n",    TimeSteps_Between_Prints);
  printf(         "Print_Forces =                %u\n",    Print_Forces);
  printf(         "Parallel execution =          ");
  #if defined(_OPENMP)
    printf(       "1\n");
    printf(       "Number of procs =             %u\n",omp_get_num_procs());
  #else
    printf(       "0\n");
  #endif

  // Are we running a new simulation or loading an existing one?
  if(Load_Data_From_File == 1) {
    timer1 = clock();

    // If loading an existing simulation, read in Particle arrays from file
    printf(       "\nLoading from file....");
    Data_Dump::Load_Simulation(&Arrays, Num_Arrays);

    // Now set up the time step counters
    Time_Step_Index = new unsigned[Num_Arrays];
    for(i = 0; i < Num_Arrays; i++)
      Time_Step_Index[i] = 0;

    timer1 = clock() - timer1;

    #if defined(_OPENMP)
      printf(     "Done!\ntook %lf s\n", timer1);
    #else
      unsigned long MS_Load = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);

      printf(       "Done!\ntook %lu ms\n", MS_Load);
    #endif
  } //   if(Load_Data_From_File == 1) {

  else if(Load_Data_From_File == 0) {
    // Use arrays defined in Simulation.h
    Use_Arrays_From_Code();

    // First, allocate the array of Bodys, time step counters
    Arrays = new Body[Num_Arrays];
    Time_Step_Index = new unsigned[Num_Arrays];
    for(i = 0; i < Num_Arrays; i++)
      Time_Step_Index[i] = 0;

    // Now set up each array using the paramaters in Simulation.h
    for(m = 0; m < Num_Arrays; m++) {
      // Set Partilce Array's name
      Arrays[m].Set_Name(Names[m]);

      // Set inter particle spacing
      Arrays[m].Set_Inter_Particle_Spacing(IPS[m]);

      // Now set the ith Body's material
      Arrays[m].Set_Material(Materials[m]);

      // Now set wheather or not the ith Body is damagable
      Arrays[m].Set_Damageable(Is_Damagable[m]);

      // Now set other Body members.
      Set_Body_Members(Arrays[m]);

      //////////////////////////////////////////////////////////////////////////
      // Check for bad inputs!

      // A body can't both be a cuboid and be from an FEB file.
      if(Is_Cuboid[m] == true && From_FEB_File[m] == true) {
        printf("A body can't be read from a FEB file and designated as a cuboid... aborting\n");
        return;
      } // if(Is_Cuboid[i] == true && From_FEB_File[i] == true) {

      // A body must either be a cuboid or be from file. If it's neither, then
      // we have no way of setting it up.
      if(Is_Cuboid[m] == false && From_FEB_File[m] == false) {
        printf("Error! All bodies must be from a FEB file or a cuboid.  Aborting\n");
        return;
      } // if(Is_Cuboid[m] == false && Is_Boundary == false) {

      //////////////////////////////////////////////////////////////////////////
      // Now set up the Body's particles

      /* If the body is a boundary, we need to designate it as such.
      note: this applies if the body is a cuboid, or from FEB file... anything
      can be a boundary! */
      if(Is_Boundary[m] == true)
        Arrays[m].Set_Boundary(true);

      // Set up body as a cuboid if it is a cuboid
      if(Is_Cuboid[m] == true) {
        Arrays[m].Set_Cuboid_Dimensions(Dimensions[m]);
        Setup_Cuboid(Arrays[m], m);
      } // if(Is_Cuboid[m] == true) {

      // if the body is from file, read it in
      else if(From_FEB_File[m] == true)
        Setup_FEB_Body(Arrays[m], m);

    } // for(m = 0; m < Num_Arrays; m++) {
  } // else if(Load_Data_From_File == 0) {

  // Now that the particle array is loaded, print paramaters.
  printf(         "\nRuinning with the following paramaters....\n");
  for(m = 0; m < Num_Arrays; m++)
    Arrays[m].Print_Parameters();

  //////////////////////////////////////////////////////////////////////////////
  // Run time steps

  printf(         "\nRunning %d time steps....\n",Num_Steps);

  // Cycle through time steps.
  timer1 = clock();
  printf(         "0 time steps complete\n");

  // If we are starting a new simulation (not reading one from file) then print
  // initial configuration
  if(Load_Data_From_File == false) {
    for(m = 0; m < Num_Arrays; m++)
      VTK_File::Export_Particle_Positions(Arrays[m]);

    if(Print_Forces == true)
      for(m = 0; m < Num_Arrays; m++)
        Arrays[m].Print_Particle_Forces();
  } // if(Load_Data_From_File == false) {

  // time step loop.
  #pragma omp parallel default(shared) private(i, j, k, l, m) firstprivate(Num_Arrays, Num_Steps, dt, TimeSteps_Between_Prints)
  {
  for(l = 0; l < Num_Steps; l++) {
    ////////////////////////////////////////////////////////////////////////////
    // Apply Boundary conditions
    // Note: we only apply BC's to 0th array.

    #pragma omp single nowait
      timer2 = clock();

    #pragma omp single
    for(m = 0; m < Num_Arrays; m++) {
      if(m == 0 && Time_Step_Index[0] == 0) {
        ////////////////////////////////////////////////////////////////////////
        /* Boundary conditions
        Here we set the Bc's for the six sides of the cube. The faces are named
        'Front', 'Back', 'Top', 'Bottom', 'Left' and 'Right'. We give the faces
        these names based on the following coordinate axis layout:

                                  Y
                                  |        X
                                  |      /
                                  |    /
                                  |  /
                              _ _ |/_ _ _ _ _ _ _ Z
                                 /|
                               /  |

        The Normal vector to the 'Front' and 'Back' faces point in the +X and -X
        directions respectivly. The Normal vector to the 'Top' and 'Bottom' faces
        point in the +Y and -Y directions respectivly. Finally, the Normal vector
        to the 'Left' and 'Right' faces point in the -Z and +Z directions
        respectivly. */

        // Establish side lengths (we assume Array[m] is a cuboid)
        unsigned X_SIDE_LENGTH = Arrays[m].Get_X_SIDE_LENGTH();
        unsigned Y_SIDE_LENGTH = Arrays[m].Get_Y_SIDE_LENGTH();
        unsigned Z_SIDE_LENGTH = Arrays[m].Get_Z_SIDE_LENGTH();

        // Front face (i = 0)
        i = 0;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
          }}

        // Back face (i = X_SIDE_LENGTH-1)
        i = X_SIDE_LENGTH-1;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
          }}

        // Bottom face (j = 0)
        j = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[1] = 0;
            //(Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V = {0,-30,0};
          }}

        // Top face (j = y_Side_len-1)
        j = Y_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            //(Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V = {0,30,0};
          }}

        // Left face (k = 0)
        k = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
          }}

        // Right face (k = Z_SIDE_LENGTH-1)
        k = Z_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
          }}

      } // if(m == 0)

      // Needle BC's
      else if(m == 1) {
        unsigned Array_m_Num_Particles = (Arrays[m]).Get_Num_Particles();
        for(i = 0; i < Array_m_Num_Particles; i++)
          if((Arrays[m])[i].Get_X()[1] > 17.)    // if y component is above threshold, press it.
            (Arrays[m])[i].V = {0, -50, 0};
      } // else if(m == 1) {

    } // for(m = 0; m < Num_Arrays; m++)

    #pragma omp single nowait
      update_BC_timer += clock() - timer2;



    ////////////////////////////////////////////////////////////////////////////
    // Update Stress tensor (P)

    #pragma omp single nowait
      timer2 = clock();

    for(m = 0; m < Num_Arrays; m++) {
      // Note: We don't update P for Bodys that are boundaries
      if(Arrays[m].Get_Boundary() == true)
        continue;
      else
        /* Update each Particles's P tensor.
        We only update P when the mth Body's counter is zero. Note
        That the Update_P method has an orphaned for loop (and takes care of
        removing damaged particles in parallel, damaged particles are not
        removed until every particle's P tensor has been updated. This makes
        the code parallelizable and determinstic) */
        if(Time_Step_Index[m] == 0)
          Arrays[m].Update_P(Steps_Per_Update[m]*dt);
    } // for(m = 0; m < Num_Arrays; m++) {

    #pragma omp single nowait
      update_P_timer += clock() - timer2;



    ////////////////////////////////////////////////////////////////////////////
    // Contact
    /* Here we enable particle-particle contact. To do this, we cycle through
    each Body. For the mth array, we check if any of its particles are
    in contact with any of the partilces in the ith array for i > m. We only
    use i > m so that we only run the contact algorythm on each part of
    Bodys once. Further, we only calculate the contact forces for the
    mth Body if that partilce_array is being updated this time step. */

    #pragma omp single nowait
      timer2 = clock();

    /* First, we need to set each particle's contact force to zero. It should be
    noted that we only do this for a particular Body if that array
    is updating it's position this turn. Otherwise, since the force won't be
    used for anything, there's no reason to waste CPU cycles setting that
    array's particle's contact forces to zero. */
    for(m = 0; m < Num_Arrays; m++) {
      unsigned Num_Particles = (Arrays[m]).Get_Num_Particles();

      if(Time_Step_Index[m] == 0)
        #pragma omp for
        for(i = 0; i < Num_Particles; i++) {
          (Arrays[m])[i].Force_Contact = {0,0,0};
          (Arrays[m])[i].Force_Friction = {0,0,0};
        } // for(i = 0; i < (Arrays[m]).Get_Num_Particles(); i++) {
    } // for(m = 0; m < Num_Arrays; m++) {

    /* Now we can apply the contact algorythm. Note that this must be applied
    every time step no matter what (so that bodies that update each step can
    are proprly updated/have the right forces applied each timestpe) */
    for(m = 0; m < Num_Arrays; m++)
      for(i = m + 1; i < Num_Arrays; i++)
        Body::Contact(Arrays[i], Arrays[m]);

    #pragma omp single nowait
      contact_timer += clock() - timer2;



    ////////////////////////////////////////////////////////////////////////////
    // Update Position (x)

    #pragma omp single nowait
      timer2 = clock();

    for(m = 0; m < Num_Arrays; m++) {
      // Note: we don't update P for Bodys that are boundaries
      if(Arrays[m].Get_Boundary() == true)
        continue;
      else {
        /* We only want to update x (the traditional way) if we're on a timestep
        where the mth Body gets updated. Suppose that the mth particle
        array only updates once every k steps (meaning that Stpes_Between_Update[m] = k)
        on the 0th step, the mth particle arrays's counter is zero. After
        each step it increments. On the kth step, its counter reaches k and the
        counter gets truncaed back to zero. Therefore, every k steps the mth
        Body's counter will be zero. Thus, we use a 0 counter
        as an indicator that we should update this Body. */
        if(Time_Step_Index[m] == 0) {
          /* First, update the 'F_Index' for the current Body. This
          controls which member of each particle's 'F' array is the 'newest'. */
          #pragma omp single
            Arrays[m].Increment_F_Index();

          // Now update the position of each particle in this body.
          Arrays[m].Update_x(dt);
        } //         if(Time_Step_Index[m] == 0) {
        else {
          /* If we're not on an update step, then we'll let this body continue
          accelerating at whatever acceleration it attained after the last
          time step. */
          unsigned Num_Particles = (Arrays[m]).Get_Num_Particles();

          #pragma omp for
          for(i = 0; i < Num_Particles; i++) {
            if((Arrays[m])[i].Get_D() >= 1)
              continue;

            (Arrays[m])[i].x += dt*(Arrays[m])[i].V;       // x_i+1 = x_i + dt*v_(i+1/2)           : mm Vector
            (Arrays[m])[i].V += dt*(Arrays[m])[i].a;      // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s Vector
          } // for(i = 0; i < (Arrays[m]).Get_Num_Particles(); i++) {
        } // else {
      } // else {
    } // for(m = 0; m < Num_Arrays; m++) {

    #pragma omp single nowait
      update_x_timer += clock() - timer2;



    ////////////////////////////////////////////////////////////////////////////
    // Update each time step counter
    /* Here we increment each Body's counter. If a particular counter
    reaches its limit (the value of Steps_Per_Update[m]) then we set that
    counter to zero (reset the counter). */

    #pragma omp single
    for(m = 0; m < Num_Arrays; m++) {
      Time_Step_Index[m]++;

      if(Time_Step_Index[m] == Steps_Per_Update[m])
        Time_Step_Index[m] = 0;
    } // for(m = 0; m < Num_Arrays; m++) {



    ////////////////////////////////////////////////////////////////////////////
    // Print to file
    #pragma omp single nowait
      timer2 = clock();

    if((l+1)%TimeSteps_Between_Prints == 0) {
      #pragma omp single nowait
        printf(     "%d time steps complete\n",l+1);

      #pragma omp for nowait
      for(m = 0; m < Num_Arrays; m++ )
        VTK_File::Export_Particle_Positions(Arrays[m]);

      if(Print_Forces == true) {
        #pragma omp for nowait
        for(m = 0; m < Num_Arrays; m++ )
          Arrays[m].Print_Particle_Forces();
      } // if(Print_Forces == true) {

      if(Print_Net_Force == true) {
        #pragma omp single nowait
        Arrays[1].Print_Net_External_Force(l+1);
      } // if(Print_Net_Force == true) {
    } // if((k+1)%100 == 0) {

    #pragma omp single nowait
      Print_timer += clock()-timer2;
  } // for(l = 0; l < Num_Steps; l++) {
  } // #pragma omp parallel
  printf(         "Done!\n");
  timer1 = clock()-timer1;

  // If saving is enabled, Dump particle data to file
  if(Save_Data_To_File == 1)
    Data_Dump::Save_Simulation(Arrays, Num_Arrays);

  // Print timing data
  #if defined(_OPENMP)
    printf(       "\nIt took %lf s to perform %u Particle time steps \n",timer1, Num_Steps);
    printf(       "%lf s to update BC's\n", update_BC_timer);
    printf(       "%lf s to update P\n", update_P_timer);
    printf(       "%lf s to update Contact\n", contact_timer);
    printf(       "%lf s to update x\n", update_x_timer);
    printf(       "%lf s to print data to files\n", Print_timer);

  #else
    unsigned long MS_Iter,                                           // These are used to store the number of miliseconds that
                  MS_BC,                                             // it took to execute each of the major operations in
                  MS_P,                                              // the code. These are only used the the code is executed
                  MS_Contact,                                        // sequentially
                  MS_x,
                  MS_Print;

    MS_Iter = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
    MS_BC = (unsigned long)((double)update_BC_timer / (double)CLOCKS_PER_MS);
    MS_P = (unsigned long)((double)update_P_timer / (double)CLOCKS_PER_MS);
    MS_Contact = (unsigned long)((double)contact_timer / (double)CLOCKS_PER_MS);
    MS_x = (unsigned long)((double)update_x_timer / (double)CLOCKS_PER_MS);
    MS_Print = (unsigned long)((double)Print_timer / (double)CLOCKS_PER_MS);

    printf(         "\nIt took %lu ms to perform %u Particle time steps \n",MS_Iter, Num_Steps);
    printf(         "%lu ms to update BC's\n", MS_BC);
    printf(         "%lu ms to update P\n", MS_P);
    printf(         "%lu ms to update Contact\n", MS_Contact);
    printf(         "%lu ms to update x\n", MS_x);
    printf(         "%lu ms to print data to files\n", MS_Print);
  #endif

  delete [] Arrays;
} // void Simulation(void) {



void Simulation::Setup_Cuboid(Body & Particles, const unsigned m) {
  unsigned i,j,k;
  clock_t timer1;

  // Particle paramaters
  const double IPS = Particles.Get_Inter_Particle_Spacing();                   //        : mm
  const double Particle_Volume = IPS*IPS*IPS;                                  //        : mm^3
  const double Particle_Radius = IPS*.578;                                     //        : mm
  const double Particle_Mass = Particle_Volume*Particles.Get_density();        //        : g

  // Furst, let's get number of partilces in the Bodys
  const unsigned X_SIDE_LENGTH = Particles.Get_X_SIDE_LENGTH();
  const unsigned Y_SIDE_LENGTH = Particles.Get_Y_SIDE_LENGTH();
  const unsigned Z_SIDE_LENGTH = Particles.Get_Z_SIDE_LENGTH();

  // Vectors to hold onto Parameters
  Vector X, x;
  Vector V = Initial_Velocity[m];                          // Initial_Velocity set in Simulation.h           : mm/s

  //////////////////////////////////////////////////////////////////////////////
  // Set up particles
  printf(         "\nGenerating particles for %s...",Particles.Get_Name().c_str());
  timer1 = clock();

  // Set up Particles
  /* Store particles in 'Vertical Column' major 'Row' semi-major order
  A vertical column is a set of particles with the same x and z coordinates,
  while a row is a set of particles with the same y and x coordinates. This
  ordering method places particles with the same */
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
        X = {i*IPS, j*IPS, k*IPS};
        X += Offset[m];
        x = X;                                                                 //        : mm

        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Mass(Particle_Mass);  //        : g
        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Vol(Particle_Volume); //        : mm^3
        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Radius(Particle_Radius);   //   : mm
        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_X(X);                 //        : mm
        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_x(x);                 //        : mm
        Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_V(V);                 //        : mm/s
      } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  timer1 = clock()-timer1;
  #if defined(_OPENMP)
    printf(        "Done!\ntook %lf s\n",timer1);
  #else
    unsigned long MS_Gen = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
    printf(        "Done!\ntook %lums\n",MS_Gen);
  #endif

  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors (if the body is not a boundary)

  if(Particles.Get_Boundary() == false) {
    printf(         "Generating %s's neighbor lists...", Particles.Get_Name().c_str());
    timer1 = clock();
    Particles.Find_Neighbors_Cuboid();

    timer1 = clock() - timer1;
    #if defined(_OPENMP)
      printf(        "Done!\ntook %lf s\n",timer1);
    #else
      unsigned long MS_Neighbor = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
      printf(       "Done!\ntook %lums\n",MS_Neighbor);
    #endif
  } //   if(Particles.Get_Boundary() == false) {

  /*
  // Damage the 'cut'
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Particles.Remove_Damaged_Particle(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2));
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {
  */
} // void Simulation::Setup_Cuboid(Body & Particles, const unsigned m) {



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
