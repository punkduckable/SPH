#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "IO/VTK_File.h"
#include "IO/Data_Dump.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif


void Simulation::Run_Simulation(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned i,j,k,l,m;

  // Computation time measurement variables
  #if defined(_OPENMP)
    double timer1,
           timer2,
           update_BC_timer = 0,
           update_P_timer = 0,
           contact_timer = 0,
           update_x_timer = 0,
           Print_timer = 0;
  #else
    clock_t timer1,
            timer2,
            update_BC_timer = 0,
            update_P_timer = 0,
            contact_timer = 0,
            update_x_timer = 0,
            Print_timer = 0;
  #endif

  // Set up Bodys.
  Body * Bodies;                                 // Will point to the Bodies's for this simulation
  unsigned * Time_Step_Index;                    // Time step counters for each body



  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.

  //Display that the simulation has begun
  printf(         "\nRunning a Simulation...\n");
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
    #if defined(_OPENMP)
      timer1 = omp_get_wtime();
    #else
      timer1 = clock();
    #endif

    // If loading an existing simulation, read in bodies from file
    printf(       "\nLoading from file....");
    Data_Dump::Load_Simulation(&Bodies, Num_Bodies);

    // Now set up the time step counters
    Time_Step_Index = new unsigned[Num_Bodies];
    for(i = 0; i < Num_Bodies; i++) {
      Time_Step_Index[i] = 0;
    } // for(i = 0; i < Num_Bodies; i++) {


    #if defined(_OPENMP)
      timer1 = omp_get_wtime() - timer1;
      printf(     "Done!\ntook %lf s\n", timer1);
    #else
      timer1 = clock() - timer1;
      unsigned long MS_Load = (unsigned long)((double)timer1 / (double)CLOCKS_PER_MS);
      printf(       "Done!\ntook %lu ms\n", MS_Load);
    #endif
  } //   if(Load_Data_From_File == 1) {

  else if(Load_Data_From_File == 0) {
    // Use Bodies defined in Simulation.h
    Body_Needle_Set_Up();

    // First, allocate the array of Bodys, time step counters
    Bodies = new Body[Num_Bodies];
    Time_Step_Index = new unsigned[Num_Bodies];
    for(i = 0; i < Num_Bodies; i++) {  Time_Step_Index[i] = 0; }

    // Now set up each body using the paramaters in Simulation.h
    for(m = 0; m < Num_Bodies; m++) {
      // Set Partilce Body's name
      Bodies[m].Set_Name(Names[m]);

      // Set inter particle spacing
      Bodies[m].Set_Inter_Particle_Spacing(IPS[m]);

      // Now set the ith Body's material
      Bodies[m].Set_Material(Simulation_Materials[m]);

      // Now set wheather or not the ith Body is damagable
      Bodies[m].Set_Damageable(Is_Damagable[m]);

      // Now set other Body members.
      Set_Body_Members(Bodies[m]);

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
      if(Is_Boundary[m] == true) {
        Bodies[m].Set_Boundary(true);
      } // if(Is_Boundary[m] == true) {

      // Set up body as a cuboid if it is a cuboid
      if(Is_Cuboid[m] == true) {
        Bodies[m].Set_Cuboid_Dimensions(Dimensions[m]);
        Setup_Cuboid(Bodies[m], m);
      } // if(Is_Cuboid[m] == true) {

      // if the body is from file, read it in
      else if(From_FEB_File[m] == true) {
        Setup_FEB_Body(Bodies[m], m);
      } // else if(From_FEB_File[m] == true) {

    } // for(m = 0; m < Num_Bodies; m++) {
  } // else if(Load_Data_From_File == 0) {

  // Now that the body is loaded, print paramaters.
  printf(         "\nRuinning with the following paramaters....\n");
  for(m = 0; m < Num_Bodies; m++) {
    Bodies[m].Print_Parameters();
  } // for(m = 0; m < Num_Bodies; m++) {



  //////////////////////////////////////////////////////////////////////////////
  // Run time steps

  printf(         "\nRunning %d time steps....\n",Num_Steps);

  // Cycle through time steps.
  #if defined(_OPENMP)
    timer1 = omp_get_wtime();
  #else
    timer1 = clock();
  #endif

  printf(         "0 time steps complete\n");

  // If we are starting a new simulation (not reading one from file) then print
  // initial configuration
  if(Load_Data_From_File == false) {
    for(m = 0; m < Num_Bodies; m++) {
      VTK_File::Export_Particle_Positions(Bodies[m]);
    } // for(m = 0; m < Num_Bodies; m++) {

    if(Print_Forces == true) {
      for(m = 0; m < Num_Bodies; m++) {
        Bodies[m].Print_Particle_Forces();
      } // for(m = 0; m < Num_Bodies; m++) {
    } // if(Print_Forces == true) {
  } // if(Load_Data_From_File == false) {

  // time step loop.
  #pragma omp parallel default(shared) private(i, j, k, l, m) firstprivate(Num_Bodies, Num_Steps, dt, TimeSteps_Between_Prints)
  {
  for(l = 0; l < Num_Steps; l++) {

    ////////////////////////////////////////////////////////////////////////////
    // Apply Boundary conditions
    // Note: we only apply BC's to 0th body.

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        timer2 = omp_get_wtime();
      #else
        timer2 = clock();
      #endif
    } // #pragma omp single nowait

    #pragma omp single
    {
    for(m = 0; m < Num_Bodies; m++) {
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

        // Establish side lengths (we assume Bodies[m] is a cuboid)
        unsigned X_SIDE_LENGTH = Bodies[m].Get_X_SIDE_LENGTH();
        unsigned Y_SIDE_LENGTH = Bodies[m].Get_Y_SIDE_LENGTH();
        unsigned Z_SIDE_LENGTH = Bodies[m].Get_Z_SIDE_LENGTH();

        // Front face (i = 0)
        i = 0;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
          }
        }

        // Back face (i = X_SIDE_LENGTH-1)
        i = X_SIDE_LENGTH-1;
        for(j = 0; j < Y_SIDE_LENGTH; j++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[0] = 0;
          }
        }

        // Bottom face (j = 0)
        j = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            (Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[1] = 0;
            //(Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V = {0,-30,0};
          }
        }

        // Top face (j = y_Side_len-1)
        j = Y_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(k = 0; k < Z_SIDE_LENGTH; k++) {
            //(Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V = {0,30,0};
          }
        }

        // Left face (k = 0)
        k = 0;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
          }
        }

        // Right face (k = Z_SIDE_LENGTH-1)
        k = Z_SIDE_LENGTH-1;
        for(i = 0; i < X_SIDE_LENGTH; i++) {
          for(j = 0; j < Y_SIDE_LENGTH; j++) {
            (Bodies[m])[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j].V[2] = 0;
          }
        }

      } // if(m == 0)

      // Needle BC's
      else if(m == 1) {
        unsigned Array_m_Num_Particles = (Bodies[m]).Get_Num_Particles();
        for(i = 0; i < Array_m_Num_Particles; i++) {
          if((Bodies[m])[i].Get_X()[1] > 17.) {    // if y component is above threshold, press it.
            (Bodies[m])[i].V = {0, -50, 0};
          } // if((Bodies[m])[i].Get_X()[1] > 17.) {
        } // for(i = 0; i < Array_m_Num_Particles; i++) {
      } // else if(m == 1) {
    } // for(m = 0; m < Num_Bodies; m++)
    } // #pragma omp single

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        update_BC_timer += omp_get_wtime() - timer2;
      #else
        update_BC_timer += clock() - timer2;
      #endif
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update Stress tensor (P)

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        timer2 = omp_get_wtime();
      #else
        timer2 = clock();
      #endif
    } // #pragma omp single nowait

    for(m = 0; m < Num_Bodies; m++) {
      // Note: We don't update P for Bodys that are boundaries
      if(Bodies[m].Get_Boundary() == true) { continue; }

      else {
        /* Update each Particles's P tensor.
        We only update P when the mth Body's counter is zero.

        Note: the Update_P method has an orphaned for loop (and takes care of
        removing damaged particles in parallel, damaged particles are not
        removed until every particle's P tensor has been updated. This makes
        the code parallelizable and determinstic) */
        if(Time_Step_Index[m] == 0) {
          Bodies[m].Update_P(Steps_Per_Update[m]*dt);
        } // if(Time_Step_Index[m] == 0) {
      } // else
    } // for(m = 0; m < Num_Bodies; m++) {

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        update_P_timer += omp_get_wtime() - timer2;
      #else
        update_P_timer += clock() - timer2;
      #endif
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Contact
    /* Here we enable particle-particle contact. To do this, we cycle through
    each Body. For the mth body, we check if any of its particles are
    in contact with any of the partilces in the ith body for i > m. We only
    use i > m so that we only run the contact algorythm on each part of
    Bodys once. Further, we only calculate the contact forces for the
    mth Body if that body is being updated this time step. */

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        timer2 = omp_get_wtime();
      #else
        timer2 = clock();
      #endif
    } // #pragma omp single nowait

    /* First, we need to set each particle's contact force to zero. It should be
    noted that we only do this for a particular Body if that body
    is updating it's position this turn. Otherwise, since the force won't be
    used for anything, there's no reason to waste CPU cycles setting that
    bodies's particle's contact forces to zero. */
    for(m = 0; m < Num_Bodies; m++) {
      unsigned Num_Particles = (Bodies[m]).Get_Num_Particles();

      if(Time_Step_Index[m] == 0) {
        #pragma omp for
        for(i = 0; i < Num_Particles; i++) {
          (Bodies[m])[i].Force_Contact = {0,0,0};
          (Bodies[m])[i].Force_Friction = {0,0,0};
        } // for(i = 0; i < (Bodies[m]).Get_Num_Particles(); i++) {
      } // if(Time_Step_Index[m] == 0) {
    } // for(m = 0; m < Num_Bodies; m++) {

    /* Now we can apply the contact algorythm. Note that this must be applied
    every time step no matter what (so that bodies that update each step can
    are proprly updated/have the right forces applied each timestpe) */
    for(m = 0; m < Num_Bodies; m++) {
      for(i = m + 1; i < Num_Bodies; i++) {
        Body::Contact(Bodies[i], Bodies[m]);
      } // for(i = m + 1; i < Num_Bodies; i++) {
    } // for(m = 0; m < Num_Bodies; m++) {

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        contact_timer = omp_get_wtime() - timer2;
      #else
        contact_timer = clock() - timer2;
      #endif
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update Position (x)

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        timer2 = omp_get_wtime();
      #else
        timer2 = clock();
      #endif
    } // #pragma omp single nowait

    for(m = 0; m < Num_Bodies; m++) {
      // Note: we don't update P for Bodys that are boundaries
      if(Bodies[m].Get_Boundary() == true) { continue; }

      else {
        /* We only want to update x (the traditional way) if we're on a timestep
        where the mth Body gets updated. Suppose that the mth body
        only updates once every k steps (meaning that Stpes_Between_Update[m] = k)
        on the 0th step, the mth Bodies's counter is zero. After
        each step it increments. On the kth step, its counter reaches k and the
        counter gets truncaed back to zero. Therefore, every k steps the mth
        Body's counter will be zero. Thus, we use a 0 counter
        as an indicator that we should update this Body. */
        if(Time_Step_Index[m] == 0) {
          /* First, update the 'F_Index' for the current Body. This
          controls which member of each particle's 'F' array is the 'newest'. */
          #pragma omp single
            Bodies[m].Increment_F_Index();

          // Now update the position of each particle in this body.
          Bodies[m].Update_x(dt);
        } //         if(Time_Step_Index[m] == 0) {
        else {
          /* If we're not on an update step, then we'll let this body continue
          accelerating at whatever acceleration it attained after the last
          time step. */
          unsigned Num_Particles = (Bodies[m]).Get_Num_Particles();

          #pragma omp for
          for(i = 0; i < Num_Particles; i++) {
            if((Bodies[m])[i].Get_D() >= 1) { continue; }

            (Bodies[m])[i].x += dt*(Bodies[m])[i].V;       // x_i+1 = x_i + dt*v_(i+1/2)           : mm Vector
            (Bodies[m])[i].V += dt*(Bodies[m])[i].a;       // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s Vector
          } // for(i = 0; i < (Bodies[m]).Get_Num_Particles(); i++) {
        } // else {
      } // else {
    } // for(m = 0; m < Num_Bodies; m++) {

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        update_x_timer += omp_get_wtime() - timer2;
      #else
        update_x_timer += clock() - timer2;
      #endif
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update each time step counter
    /* Here we increment each Body's counter. If a particular counter
    reaches its limit (the value of Steps_Per_Update[m]) then we set that
    counter to zero (reset the counter). */

    #pragma omp single
    for(m = 0; m < Num_Bodies; m++) {
      Time_Step_Index[m]++;

      if(Time_Step_Index[m] == Steps_Per_Update[m])
        Time_Step_Index[m] = 0;
    } // for(m = 0; m < Num_Bodies; m++) {



    ////////////////////////////////////////////////////////////////////////////
    // Print to file
    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        timer2 = omp_get_wtime();
      #else
        timer2 = clock();
      #endif
    } // #pragma omp single nowait

    if((l+1)%TimeSteps_Between_Prints == 0) {
      #pragma omp single nowait
        printf(     "%d time steps complete\n",l+1);

      #pragma omp for nowait
      for(m = 0; m < Num_Bodies; m++ ) {
        VTK_File::Export_Particle_Positions(Bodies[m]);
      } // for(m = 0; m < Num_Bodies; m++ ) {

      if(Print_Forces == true) {
        #pragma omp for nowait
        for(m = 0; m < Num_Bodies; m++ ) {
          Bodies[m].Print_Particle_Forces();
        } // for(m = 0; m < Num_Bodies; m++ ) {
      } // if(Print_Forces == true) {

      if(Print_Net_Force == true) {
        #pragma omp single nowait
        Bodies[1].Print_Net_External_Force(l+1);
      } // if(Print_Net_Force == true) {
    } // if((k+1)%100 == 0) {

    #pragma omp single nowait
    {
      #if defined(_OPENMP)
        Print_timer += omp_get_wtime() - timer2;
      #else
        Print_timer += clock()-timer2;
      #endif
    } // #pragma omp single nowait
  } // for(l = 0; l < Num_Steps; l++) {
  } // #pragma omp parallel
  printf(         "Done!\n");
  #if defined(_OPENMP)
    timer1 = omp_get_wtime() - timer1;
  #else
    timer1 = clock()-timer1;
  #endif

  // If saving is enabled, Dump particle data to file
  if(Save_Data_To_File == 1) {
    Data_Dump::Save_Simulation(Bodies, Num_Bodies);
  } // if(Save_Data_To_File == 1) {

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

  delete [] Bodies;
} // void Simulation(void) {



void Startup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {
  // Fill me in!
} // void Startup_Simulation(Body ** Bodies, unsigned ** Time_Step_Index) {
