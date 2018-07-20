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
  //        contact_timer = 0,
          update_x_timer = 0,
          Print_timer = 0;
  unsigned long MS_Iter,
                MS_BC,
                MS_P,
  //              MS_Contact,
                MS_x,
                MS_Print;                                            // Timers (store number of MS for each operation)

  // Partile bodies
  Particle *Body;
  //Particle *Boundary;

  unsigned int Num_Particles_Body;
  //const unsigned int Num_Particles_Boundary = 4*X_SIDE_LENGTH*Z_SIDE_LENGTH;

  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.

  //Display that the simulation has begun
  printf(         "\nRunning a simulation...\n");
  printf(         "Load_Data_From_File =         %u\n",Load_Data_From_File);
  printf(         "Save_Data_To_File =           %u\n",Save_Data_To_File);

  // Are we running a new simulation or loading an existing one?
  if(Load_Data_From_File == 1) {
    unsigned long MS_Load;
    timer1 = clock();

    // If loading an existing simulation, read in 'Particle_File'
    printf(       "Loading from file....\n");
    Body = Data_Dump::Load_Data_From_File(Num_Particles_Body);

    timer1 = clock() - timer1;
    MS_Load = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
    printf(       "Done! Took %lu ms\n", MS_Load);
  } //   if(Load_Data_From_File == 1) {
  else if(Load_Data_From_File == 0) {
    // Set up default static values (for the particle class)
    Set_Static_Particle_Members();

    /* With the static members set, let's create a particle's array.*/
    Num_Particles_Body = X_SIDE_LENGTH*Y_SIDE_LENGTH*Z_SIDE_LENGTH;
    Body = new Particle[Num_Particles_Body];
    //Boundary = new Particle[Num_Particles_Boundary];

    // Now let's setup the particle's array (this function sets each particle's
    // position, ID, etc... It also finds and sets each particle's neighbors)
    Set_Up_Body(Body, Num_Particles_Body, Particle::Inter_Particle_Spacing);
    //Set_Up_Boundary(Boundary, Num_Particles_Boundary, Particle::Inter_Particle_Spacing);
  } // else if(Load_Data_From_File == 0) {

  // Now that the particle array is loaded, print paramaters.
  printf(         "Ruinning with the following paramaters....\n");
  printf(         "X side length:                %u\n", X_SIDE_LENGTH);
  printf(         "Y side length:                %u\n", Y_SIDE_LENGTH);
  printf(         "Z side length:                %u\n", Z_SIDE_LENGTH);
  printf(         "Inter particle spacing:       %lf\n", Particle::Inter_Particle_Spacing);
  printf(         "h:                            %lf\n", Particle::h);
  printf(         "Support Radius:               %u\n", Particle::Support_Radius);
  printf(         "Shape Function Amplitude:     %lf\n", Particle::Shape_Function_Amp);
  printf(         "Lame:                         %lf\n", Particle::Lame);
  printf(         "mu0 (Shear modulus):          %lf\n", Particle::mu0);
  printf(         "mu (Viscosity):               %lf\n", Particle::mu);
  printf(         "E (Young's modulus):          %lf\n", Particle::E);
  printf(         "Tau (Damage rate):            %lf\n\n", Particle::Tau);

  //////////////////////////////////////////////////////////////////////////////
  // Run time steps
  printf(         "Running particle time steps....\n");

  // Print initial data
  timer1 = clock();
  printf(         "0 time steps complete\n");
  VTK_File::Export_Pariticle_Positions(Num_Particles_Body, Body);
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
      }
    }

    // Back face (i = X_SIDE_LENGTH-1)
    i = X_SIDE_LENGTH-1;
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
      }
    }

    // Bottom face (j = 0)
    j = 0;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].Set_V({0,-20,0});
      }
    }

    // Top face (j = y_Side_len-1)
    j = Y_SIDE_LENGTH-1;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Body[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].Set_V({0,20,0});
      }
    }

    // Left face (k = 0)
    k = 0;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
      }
    }

    // Right face (k = Z_SIDE_LENGTH-1)
    k = Z_SIDE_LENGTH-1;
    for(i = 0; i < X_SIDE_LENGTH; i++) {
      for(j = 0; j < Y_SIDE_LENGTH; j++) {
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
    update_P_timer += clock() - timer2;

    ////////////////////////////////////////////////////////////////////////////
    // Detect contact
    //timer2 = clock();
    //Contact(Body, Num_Particles_Body, Boundary, Num_Particles_Boundary, h);
    //contact_timer += clock() - timer2;

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
    update_x_timer += clock() - timer2;

    // Print to file evert 100th iteration
    timer2 = clock();
    if((l+1)%100 == 0) {
      printf(     "%d time steps complete\n",l+1);
      VTK_File::Export_Pariticle_Positions(Num_Particles_Body, Body);

      //Particle_Debugger::Export_Pariticle_Forces(Num_Particles_Body, Body);
    } // if((k+1)%100 == 0) {
    Print_timer += clock()-timer2;
  } // for(l = 0; l < Num_Steps; l++) {
  printf(         "Done!\n\n");
  timer1 = clock()-timer1;

  // If saving is enabled, Dump particle data to file
  if(Save_Data_To_File == 1)
    Data_Dump::Print_Data_To_File(Body, Num_Particles_Body);

  // Print timing data
  MS_Iter = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
  MS_BC = (unsigned long)((double)update_BC_timer / (double)CLOCKS_PER_MS);
  MS_P = (unsigned long)((double)update_P_timer / (double)CLOCKS_PER_MS);
  //MS_Contact = (unsigned long)((double)contact_timer / (double)CLOCKS_PER_MS);
  MS_x = (unsigned long)((double)update_x_timer / (double)CLOCKS_PER_MS);
  MS_Print = (unsigned long)((double)Print_timer / (double)CLOCKS_PER_MS);

  printf(         "It took %lu ms to perform %u Particle time steps \n",MS_Iter, Num_Steps);
  printf(         "%lu ms to update BC's\n", MS_BC);
  printf(         "%lu ms to update P\n", MS_P);
  //printf(         "%lu ms to update Contact\n", MS_Contact);
  printf(         "%lu ms to update x\n", MS_x);
  printf(         "%lu ms to print data to files\n", MS_Print);

  // Free memory
  delete [] Body;
  //delete [] Boundary;

} // void Simulation(void) {

void Simulation::Set_Up_Body(Particle * Body, const unsigned int Num_Particles_Body, const double IPS) {
  unsigned long MS_Gen,
                MS_Neighbor;
  clock_t timer1;
  unsigned int i,j,k;

  // Set particle paramaters
  double Particle_Volume = IPS*IPS*IPS;                                      //        : mm^3
  double Particle_Radius = IPS*.578;                                         //        : mm
  double Particle_Density = 1;                                               //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                   //        : g

  Vector X, x, V;

  //////////////////////////////////////////////////////////////////////////////
  // Set up particles
  printf(         "Generating particles....\n");
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

        X *= IPS;                                                            //        : mm
        x = X;                                                               //        : mm
        V = {0.,0.,0.};                                                      //        : mm/s

        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_ijk(i,j,k);      // Remove if not using cuboid
        Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j].Set_ID(i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j);
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
  printf(         "Done! took %lums\n\n",MS_Gen);

  //////////////////////////////////////////////////////////////////////////////
  // Set up Neighbors

  printf(         "Generating Neighbor lists....\n");
  timer1 = clock();
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      for(k = 0; k < Z_SIDE_LENGTH; k++) {
        Particle_Helpers::Find_Neighbors_Box(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + j], Body, X_SIDE_LENGTH, Y_SIDE_LENGTH, Z_SIDE_LENGTH);
      } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
    } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  // Damage the 'cut'
  for(i = 0; i < 1; i++) {                     // Depth of cut
    for(k = 0; k < Z_SIDE_LENGTH; k++) {       // Length of cut
      Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)].Set_D(1);
      Particle_Helpers::Remove_Damaged_Particle(Body[i*(Y_SIDE_LENGTH*Z_SIDE_LENGTH) + k*Y_SIDE_LENGTH + (Y_SIDE_LENGTH/2)], Body);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < 3; i++) {

  timer1 = clock() - timer1;
  MS_Neighbor = (unsigned long)(((float)timer1)/((float)CLOCKS_PER_MS));
  printf(         "Done! took %lums\n\n",MS_Neighbor);

} // void Set_Up_Body(Particle * Body, const unsigned int Num_Particles, Body, const double IPS) {

void Simulation::Set_Up_Boundary(Particle * Boundary, const unsigned int Num_Particles_Boundary, const double IPS) {
  unsigned int i,k;

  // Set particle paramaters
  double Particle_Volume = IPS*IPS*IPS;                                        //        : mm^3
  double Particle_Radius = IPS*.578;                                           //        : mm
  double Particle_Density = 1;                                                 //        : g/mm^3
  double Particle_Mass = Particle_Volume*Particle_Density;                     //        : g

  Vector X, x, V;

  // Set up Boundary
  for(i = 0; i < 2*X_SIDE_LENGTH; i++) {
    for(k = 0; k < 2*Z_SIDE_LENGTH; k++) {
      X = {(double)i - .5*(double)X_SIDE_LENGTH, -15. , (double)k - .5*(double)Z_SIDE_LENGTH};
      X *= IPS;
      x = X;
      V = {0.,0.,0.};

      Boundary[i*2*Z_SIDE_LENGTH + k].Set_ijk(i,1,k);
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_ID(i*Z_SIDE_LENGTH + k);
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Mass(Particle_Mass);                 //        : g
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Vol(Particle_Volume);                //        : mm^3
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_Radius(Particle_Radius);             //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_X(X);                                //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_x(x);                                //        : mm
      Boundary[i*2*Z_SIDE_LENGTH + k].Set_V(V);                                //        : mm/s
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {
} // void Simulation::Set_Up_Boundary(Particle * Boundary, const unsigned int Num_Particles_Boundary, const double IPS) {

#endif
