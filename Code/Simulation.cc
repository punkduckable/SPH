#if !defined(SIMULATION_SOURCE)
#define SIMULATION_SOURCE

#include "Simulation.h"

void Simulation::Run_Simulation(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned int i,j,k,l,m;

  // Computation time measurement variables
  clock_t timer1,
          timer2,
          update_BC_timer = 0,
          update_P_timer = 0,
          contact_timer = 0,
          update_x_timer = 0,
          Print_timer = 0;
  unsigned long MS_Iter,
                MS_BC,
                MS_P,
                MS_Contact,
                MS_x,
                MS_Print;                                            // Timers (store number of MS for each operation)

  // Set up Particle_Arrays.
  Particle_Array * Arrays;

  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.

  //Display that the simulation has begun
  printf(         "\nRunning a simulation...\n");
  printf(         "Load_Data_From_File =         %u\n",    Load_Data_From_File);
  printf(         "Save_Data_To_File =           %u\n",    Save_Data_To_File);
  printf(         "TimeSteps_Between_Prints =    %u\n",    TimeSteps_Between_Prints);
  printf(         "Print_Forces =                %u\n",    Print_Forces);

  // Are we running a new simulation or loading an existing one?
  if(Load_Data_From_File == 1) {
    unsigned long MS_Load;
    timer1 = clock();

    // If loading an existing simulation, read in Particle arrays from file
    printf(       "\nLoading from file....");
    Data_Dump::Load_Simulation(&Arrays, Num_Arrays);

    timer1 = clock() - timer1;
    MS_Load = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
    printf(       "Done!\ntook %lu ms\n", MS_Load);
  } //   if(Load_Data_From_File == 1) {

  else if(Load_Data_From_File == 0) {
    // Use arrays defined in Simulation.h
    Use_Arrays_From_Code();

    // First, allocate the array of Particle_Arrays
    Arrays = new Particle_Array[Num_Arrays];

    // Now set up each array using the paramaters in Simulation.h
    for(m = 0; m < Num_Arrays; m++) {
      // Set Partilce Array's name
      Arrays[m].Set_Name(Names[m]);

      // Now set the ith Particle_Array's members.
      Set_Particle_Array_Members(Arrays[m]);

      // Now set the ith Particle_Array's material
      Arrays[m].Set_Material(Materials[m]);

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
      // Now set up the Particle_Array's particles

      // Set up body as a cuboid if it is a cuboid
      if(Is_Cuboid[m] == true) {
        Arrays[m].Set_Cuboid_Dimensions(Dimensions[m]);
        Setup_Cuboid(Arrays[m]);
      } // if(Is_Cuboid[m] == true) {

      // if the body is from file, read it in
      else if(From_FEB_File[m] == true)
        Setup_FEB_Body(Arrays[m], Names[m]);

      /* If the body is a boundary, we need to designate it as such.
      note: this applies if the body is a cuboid, or from FEB file... anything
      can be a boundary! Thus, we apply this check after the setup proess */
      if(Is_Boundary[m] == true)
        Arrays[m].Set_Boundary(true);

    } // for(m = 0; m < Num_Arrays; m++) {
  } // else if(Load_Data_From_File == 0) {

  // Now that the particle array is loaded, print paramaters.
  printf(         "\nRuinning with the following paramaters....\n");
  for(m = 0; m < Num_Arrays; m++)
    Arrays[m].Print_Parameters();

  //////////////////////////////////////////////////////////////////////////////
  // Run time steps

  unsigned int X_SIDE_LENGTH;
  unsigned int Y_SIDE_LENGTH;
  unsigned int Z_SIDE_LENGTH;
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
        Particle_Debugger::Export_Particle_Forces(Arrays[m]);
  } // if(Load_Data_From_File == false) {

  // time step loop.
  for(l = 0; l < Num_Steps; l++) {
    // Cycle through the particle arrays, apply BC's to 0th array.
    timer2 = clock();
    for(m = 0; m < Num_Arrays; m++) {
      if(m == 0) {
        ////////////////////////////////////////////////////////////////////////////
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
        X_SIDE_LENGTH = Arrays[m].Get_X_SIDE_LENGTH();
        Y_SIDE_LENGTH = Arrays[m].Get_Y_SIDE_LENGTH();
        Z_SIDE_LENGTH = Arrays[m].Get_Z_SIDE_LENGTH();

        // Front face (i = 0)
        i = 0;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0; }}

        // Back face (i = X_SIDE_LENGTH-1)
        i = X_SIDE_LENGTH-1;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0; }}

        // Bottom face (j = 0)
        j = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[1] = 0; }}

        // Top face (j = y_Side_len-1)
        j = Y_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) { }}

        // Left face (k = 0)
        k = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0; }}

        // Right face (k = Z_SIDE_LENGTH-1)
        k = Z_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Arrays[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0; }}

      } // if(m == 0)
    } // for(m = 0; m < Num_Arrays; m++)
    update_BC_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Update each Particle_Array's Stress tensors

    timer2 = clock();
    for(m = 0; m < Num_Arrays; m++) {
      // Note: We don't update P for Particle_Arrays that are boundaries
      if(Arrays[m].Get_Boundary() == true)
        continue;
      else
        for(i = 0; i < (Arrays[m]).Get_Num_Particles(); i++)
          Particle_Helpers::Update_P((Arrays[m])[i], Arrays[m], dt);
    } // for(m = 0; m < Num_Arrays; m++) {
    update_P_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Contact

    timer2 = clock();
    Particle_Helpers::Contact(Arrays[0], Arrays[1]);
    contact_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Update each Particle_Array's positions

    timer2 = clock();
    for(m = 0; m < Num_Arrays; m++) {
      // Note: we don't update P for Particle_Arrays that are boundaries
      if(Arrays[m].Get_Boundary() == true)
        continue;
      else
        for(i = 0; i < (Arrays[m]).Get_Num_Particles(); i++)
          Particle_Helpers::Update_x((Arrays[m])[i], Arrays[m], dt);
    } // for(m = 0; m < Num_Arrays; m++) {
    update_x_timer += clock() - timer2;

    // Print to file
    timer2 = clock();
    if((l+1)%TimeSteps_Between_Prints == 0) {
      printf(     "%d time steps complete\n",l+1);
      for(m = 0; m < Num_Arrays; m++ )
        VTK_File::Export_Particle_Positions(Arrays[m]);

      if(Print_Forces == true)
        for(m = 0; m < Num_Arrays; m++ )
          Particle_Debugger::Export_Particle_Forces(Arrays[m]);
    } // if((k+1)%100 == 0) {
    Print_timer += clock()-timer2;
  } // for(l = 0; l < Num_Steps; l++) {
  printf(         "Done!\n");
  timer1 = clock()-timer1;

  // If saving is enabled, Dump particle data to file
  if(Save_Data_To_File == 1)
    Data_Dump::Save_Simulation(Arrays, Num_Arrays);

  // Print timing data
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

  delete [] Arrays;
} // void Simulation(void) {



void Simulation::Setup_Cuboid(Particle_Array & Particles) {
  unsigned int i,j,k;

  unsigned long MS_Gen,
                MS_Neighbor;
  clock_t timer1;

  // Particle paramaters
  const double IPS = Particles.Get_Inter_Particle_Spacing();                   //        : mm
  const double Particle_Volume = IPS*IPS*IPS;                                  //        : mm^3
  const double Particle_Radius = IPS*.578;                                     //        : mm
  const double Particle_Density = 1;                                           //        : g/mm^3
  const double Particle_Mass = Particle_Volume*Particle_Density;               //        : g

  // Furst, let's get number of partilces in the Particle_Arrays
  const unsigned int X_SIDE_LENGTH = Particles.Get_X_SIDE_LENGTH();
  const unsigned int Y_SIDE_LENGTH = Particles.Get_Y_SIDE_LENGTH();
  const unsigned int Z_SIDE_LENGTH = Particles.Get_Z_SIDE_LENGTH();

  // Vectors to hold onto Parameters
  Vector X, x, V;
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
        x = X;                                                               //        : mm
        V = {0.,0.,0.};                                                      //        : mm/s

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
  MS_Gen = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
  printf(         "Done!\ntook %lums\n",MS_Gen);

  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors

  printf(         "Generating %s's neighbor lists...", Particles.Get_Name().c_str());
  timer1 = clock();
  for(i = 0; i < X_SIDE_LENGTH; i++)
    for(j = 0; j < Y_SIDE_LENGTH; j++)
      for(k = 0; k < Z_SIDE_LENGTH; k++)
        Particle_Helpers::Find_Neighbors_Box(Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j], Particles);

  // Damage the 'cut'
  /*
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Particle_Helpers::Remove_Damaged_Particle(Particles[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)], Particles);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {
  */

  timer1 = clock() - timer1;
  MS_Neighbor = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
  printf(         "Done!\ntook %lums\n",MS_Neighbor);
} // void Simulation::Setup_Cuboid(Particle_Array & Particles) {



void Simulation::Setup_FEB_Body(Particle_Array & FEB_Body, const std::string & File_Name) {
  // First, we need to know how many particles we have, and the reference
  // position of each of the particles.
  Vector * X = NULL;
  unsigned int Num_Particles;
  FEB_File::Read_FEB_File(File_Name, &X, Num_Particles);

  printf("\nReading in Particles for %s from FEB file...\n", FEB_Body.Get_Name().c_str());

  // Now we can set up the body
  FEB_Body.Set_Num_Particles(Num_Particles);

  // Now we can cycle through the particles, setting up each particle.
  const double IPS = FEB_Body.Get_Inter_Particle_Spacing();                    //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Density = 1;                                                 //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                     //        : g

  Vector V{0,0,0};

  for(unsigned int i = 0; i < Num_Particles; i++) {
    FEB_Body[i].Set_Mass(Particle_Mass);
    FEB_Body[i].Set_Vol(Particle_Volume);
    FEB_Body[i].Set_Radius(Particle_Radius);
    FEB_Body[i].Set_X(X[i]);
    FEB_Body[i].Set_x(X[i]);
    FEB_Body[i].Set_V(V);
  } //   for(unsigned int i = 0; i < Num_Particles; i++) {

  // Now set up neighbors.
  printf("Setting up neighbors for %s...\n",FEB_Body.Get_Name().c_str());
  Particle_Helpers::Find_Neighbors(FEB_Body);
  printf("Done!\n");
} // void Simulation::Setup_FEB_Body(Particle_Array & FEB_Body, const std::string & File_Name) {

#endif
