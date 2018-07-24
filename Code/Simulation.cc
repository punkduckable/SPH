#if !defined(SIMULATION_SOURCE)
#define SIMULATION_SOURCE

#include "Simulation.h"

void Simulation::Run_Simulation(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned int i,j,k,l;

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

  // Partile bodies
  Particle_Array Body;
  //Particle_Array Boundary;
  Particle_Array Needle;

  unsigned int Num_Particles_Body;
  //unsigned int Num_Particles_Boundary;
  unsigned int Num_Particles_Needle;

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

    // If loading an existing simulation, read in 'Particle_File'
    printf(       "\nLoading from file....");
    Data_Dump::Load_Particle_Array_From_File(Body);

    Num_Particles_Body = Body.Get_Num_Particles();

    timer1 = clock() - timer1;
    MS_Load = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
    printf(       "Done!\n");
    printf(       "took %lu ms\n", MS_Load);
  } //   if(Load_Data_From_File == 1) {
  else if(Load_Data_From_File == 0) {
    // First, let's calculate the number of particles in our two arrays.
    Num_Particles_Body = X_SIDE_LENGTH*Y_SIDE_LENGTH*Z_SIDE_LENGTH;
    //Num_Particles_Boundary = 4*X_SIDE_LENGTH*Z_SIDE_LENGTH;

    // Now set up the particle arrays
    Body.Set_Num_Particles(Num_Particles_Body);
    //Boundary.Set_Num_Particles(Num_Particles_Boundary);

    // Set up Body, Boundary.
    Set_Particle_Array_Members(Body);
    //Set_Particle_Array_Members(Boundary);

    // Now let's setup the particle's array (this function sets each particle's
    // position, ID, etc... It also finds and sets each particle's neighbors)
    SetUp_Body(Body);
    //SetUp_Boundary(Boundary);
    SetUp_FEB_Body(Needle, "Needle.feb");
  } // else if(Load_Data_From_File == 0) {

  // Now that the particle array is loaded, print paramaters.
  printf(         "\nRuinning with the following paramaters....\n");
  printf(         "X side length:                %u\n",    X_SIDE_LENGTH);
  printf(         "Y side length:                %u\n",    Y_SIDE_LENGTH);
  printf(         "Z side length:                %u\n",    Z_SIDE_LENGTH);
  printf(         "Inter particle spacing:       %lf\n",   Body.Get_Inter_Particle_Spacing());
  printf(         "h:                            %lf\n",   Body.Get_h());
  printf(         "Support Radius:               %u\n",    Body.Get_Support_Radius());
  printf(         "Shape Function Amplitude:     %lf\n",   Body.Get_Shape_Function_Amplitude());
  printf(         "Lame:                         %lf\n",   Body.Get_Lame());
  printf(         "mu0 (Shear modulus):          %lf\n",   Body.Get_mu0());
  printf(         "mu (Viscosity):               %lf\n",   Body.Get_mu());
  printf(         "E (Young's modulus):          %lf\n",   Body.Get_E());
  printf(         "Tau (Damage rate):            %lf\n",   Body.Get_Tau());

  //////////////////////////////////////////////////////////////////////////////
  // Run time steps
  printf(         "\nRunning %d time steps....\n",Num_Steps);

  // Print initial data
  timer1 = clock();
  printf(         "0 time steps complete\n");
  VTK_File::Export_Pariticle_Positions(Body);
  if(Print_Forces == true)
    Particle_Debugger::Export_Pariticle_Forces(Body);
  //Particle_Debugger::Export_Pariticle_Forces(Num_Particles_Body, Body);

  // Time step loop
  for(l = 0; l < Num_Steps; l++) {
    timer2 = clock();
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

    // Front face (i = 0)
    i = 0;
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
      }
    }

    // Back face (i = X_SIDE_LENGTH-1)
    i = X_SIDE_LENGTH-1;
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
      }
    }

    // Bottom face (j = 0)
    j = 0;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[1] = 0;
      }
    }

    // Top face (j = y_Side_len-1)
    j = Y_SIDE_LENGTH-1;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
      }
    }

    // Left face (k = 0)
    k = 0;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
      }
    }

    // Right face (k = Z_SIDE_LENGTH-1)
    k = Z_SIDE_LENGTH-1;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
      }
    }

    update_BC_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Update each particle's Stress tensor
    timer2 = clock();
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          Particle_Helpers::Update_P(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j], Body, dt);
        } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
      } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(i = 0; i < X_SIDE_LENGTH; i++) {

    for(i = 0; i < Num_Particles_Needle; i++)
      Particle_Helpers::Update_P(Needle[i], Needle, dt);

    update_P_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Detect contact
    timer2 = clock();
    Particle_Helpers::Contact(Body, Needle);
    contact_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Update each particle's position
    timer2 = clock();
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          Particle_Helpers::Update_x(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j], Body, dt);
        } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
      } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(i = 0; i < X_SIDE_LENGTH; i++) {

    for(i = 0; i < Num_Particles_Needle; i++)
      Particle_Helpers::Update_x(Needle[i], Needle, dt);

    update_x_timer += clock() - timer2;

    // Print to file evert 100th iteration
    timer2 = clock();
    if((l+1)%TimeSteps_Between_Prints == 0) {
      printf(     "%d time steps complete\n",l+1);
      VTK_File::Export_Pariticle_Positions(Body);

      if(Print_Forces == true)
        Particle_Debugger::Export_Pariticle_Forces(Body);
    } // if((k+1)%100 == 0) {
    Print_timer += clock()-timer2;
  } // for(l = 0; l < Num_Steps; l++) {
  printf(         "Done!\n");
  timer1 = clock()-timer1;

  // If saving is enabled, Dump particle data to file
  if(Save_Data_To_File == 1)
    Data_Dump::Print_Particle_Array_To_File(Body);

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
} // void Simulation(void) {

void Simulation::SetUp_Body(Particle_Array & Body) {
  unsigned int i,j,k;

  unsigned long MS_Gen,
                MS_Neighbor;
  clock_t timer1;

  // Set particle paramaters
  const double IPS = Body.Get_Inter_Particle_Spacing();                        //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Density = 1;                                                 //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                     //        : g

  Vector X, x, V;
  Vector Offset{-((double)X_SIDE_LENGTH)/2., -((double)Y_SIDE_LENGTH)/2., -((double)Z_SIDE_LENGTH)/2.};
  //////////////////////////////////////////////////////////////////////////////
  // Set up particles
  printf(         "\nGenerating particles....");
  timer1 = clock();

  // Set up Body
  /* Store particles in 'Vertical Column' major 'Row' semi-major order
  A vertical column is a set of particles with the same x and z coordinates,
  while a row is a set of particles with the same y and x coordinates. This
  ordering method places particles with the same */
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
        X = {(double)i,(double)j,(double)k};
        X += Offset;

        X *= IPS;                                                            //        : mm
        x = X;                                                               //        : mm
        V = {0.,0.,0.};                                                      //        : mm/s

        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_ijk(i,j,k);      // Remove if not using cuboid
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Mass(Particle_Mass);     //        : g
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Vol(Particle_Volume);    //        : mm^3
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_Radius(Particle_Radius); //   : mm
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_X(X);                    //        : mm
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_x(x);                    //        : mm
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_V(V);                    //        : mm/s
      } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  timer1 = clock()-timer1;
  MS_Gen = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
  printf(         "Done!\n");
  printf(         "took %lums\n",MS_Gen);

  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors

  printf(         "\nGenerating Body neighbor lists....");
  timer1 = clock();
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Particle_Helpers::Find_Neighbors_Box(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j], Body, X_SIDE_LENGTH, Y_SIDE_LENGTH, Z_SIDE_LENGTH);
      } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  // Damage the 'cut'
  /*
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Particle_Helpers::Remove_Damaged_Particle(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)], Body);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {
  */
  timer1 = clock() - timer1;
  MS_Neighbor = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
  printf(         "Done!\n");
  printf(         "took %lums\n",MS_Neighbor);

} // void Simulation::SetUp_Body(Particle_Array & Body) {

void Simulation::SetUp_Boundary(Particle_Array & Boundary) {
  unsigned int i,k;

  // Set particle paramaters
  const double IPS = Boundary.Get_Inter_Particle_Spacing();                    //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Density = 1;                                                 //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                     //        : g

  Vector X, x, V;
  double X1, X2, X3;

  // Set up Boundary
  printf(         "\nGenerating boundary....");
  for(i = 0; i < 2*X_SIDE_LENGTH; i++) {
    for(k = 0; k < 2*Z_SIDE_LENGTH; k++) {
      X1 = (double)i + 5.;//- .5*(double)X_SIDE_LENGTH;
      X3 = (double)k - .5*(double)Z_SIDE_LENGTH;
      X2 = -10.;

      X = {X1, X2, X3};
      X *= IPS;
      x = X;
      V = {0.,0.,0.};

      Boundary[i*2*Z_SIDE_LENGTH + k].Set_ijk(i,1,k);
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Mass(Particle_Mass);                 //        : g
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Vol(Particle_Volume);                //        : mm^3
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Radius(Particle_Radius);             //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_X(X);                                //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_x(x);                                //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_V(V);                                //        : mm/s
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {
  printf(         "Done!\n");
} // void Simulation::SetUp_Boundary(Particle_Array & Boundary) {

void Simulation::SetUp_FEB_Body(Particle_Array & FEB_Body, const std::string & File_Name) {
  // First, we need to know how many particles we have, and the reference
  // position of each of the particles.
  Vector * X = NULL;
  unsigned int Num_Particles;

  FEB_File::Read_FEB_File(File_Name, &X, Num_Particles);

  // Now we can set up the particle array,
  FEB_Body.Set_Num_Particles(Num_Particles);
  FEB_Body.Set_Inter_Particle_Spacing(.5);
  FEB_Body.Set_Support_Radius(4);
  FEB_Body.Set_Lame(1.125);                               // Lame parameter             : Mpa
  FEB_Body.Set_mu0(.275);                                 // Shear modulus              : Mpa
  FEB_Body.Set_mu(5e-4);                                  // Viscosity                  : Mpa*s
  FEB_Body.Set_E(0.770982);                               // Youngs modulus/Hourglass stiffness   : Mpa
  FEB_Body.Set_alpha(7.5);                                // Hg control parameter       : Unitless
  FEB_Body.Set_Tau(.15);                                  // Damage rate parameter      : unitless

  // Now we can cycle through the particles, setting up each particle.
  const double IPS = 5./8.;                                                    //        : mm
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Density = 1;                                                 //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                     //        : g

  Vector x, V{0,0,0};

  for(unsigned int i = 0; i < Num_Particles; i++) {
    x = X[i];

    FEB_Body[i].Set_Mass(Particle_Mass);
    FEB_Body[i].Set_Vol(Particle_Volume);
    FEB_Body[i].Set_Radius(Particle_Radius);
    FEB_Body[i].Set_X(X[i]);
    FEB_Body[i].Set_x(x);
    FEB_Body[i].Set_V(V);
  } //   for(unsigned int i = 0; i < Num_Particles; i++) {

  // Now set up neighbors.
  Particle_Helpers::Find_Neighbors(FEB_Body);
} // void Simulation::SetUp_FEB_Body(Particle_Array & FEB_Body, const std::string & File_Name) {

#endif
