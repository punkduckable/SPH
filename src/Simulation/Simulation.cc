#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "IO/Data_Dump.h"
#include "Errors.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif


void Simulation::Run_Simulation(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned p,b,t;

  // Simulation time measurement variables
  TIME_TYPE time1,
            time2,
            simulation_time = 0,
            update_BC_time = 0,
            update_P_time = 0,
            contact_time = 0,
            update_x_time = 0,
            Print_time = 0;

  Body * Bodies;                                 // Will point to the Bodies's for this simulation
  unsigned * Time_Step_Index;                    // Time step counters for each body


  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.

  Startup_Simulation(&Bodies, &Time_Step_Index);



  //////////////////////////////////////////////////////////////////////////////
  // Run time steps
  printf(         "\nRunning %d time steps....\n",Num_Time_Steps);
  time1 = Get_Time();

  // time step loop.
  #pragma omp parallel default(shared) private(b, p, t) firstprivate(Num_Bodies, Num_Time_Steps, dt, TimeSteps_Between_Prints)
  {
  for(t = 0; t < Num_Time_Steps; t++) {
    ////////////////////////////////////////////////////////////////////////////
    // Print to file
    #pragma omp single nowait
    {
      time2 = Get_Time();
    } // #pragma omp single nowait

    if(t%TimeSteps_Between_Prints == 0) {
      #pragma omp single nowait
        printf("%d time steps complete\n",t);

      #pragma omp for nowait
      for(b = 0; b < Num_Bodies; b++) {
        try { Bodies[b].Export_Particle_Positions(); }
        catch(Exception & Error_In) {
          printf("%s\n", Error_In.what());
          abort();
        } // catch(Exception & Error_In) {
      } // for(b = 0; b < Num_Bodies; b++ ) {

      if(Print_Particle_Forces == true) {
        #pragma omp for nowait
        for(b = 0; b < Num_Bodies; b++) {
          try { Bodies[b].Export_Particle_Forces(); }
          catch(Exception & Error_In) {
            printf("%s\n", Error_In.what());
            abort();
          } // catch(Exception & Error_In) {
        } // for(b = 0; b < Num_Bodies; b++) {
      } // if(Print_Particle_Forces == true) {

      if(Print_Net_External_Forces == true) {
        #pragma omp for nowait
        for(b = 0; b < Num_Bodies; b++) {
          try { Bodies[b].Export_Net_External_Force(t); }
          catch(Exception & Error_In) {
            printf("%s\n", Error_In.what());
            abort();
          } // catch(Exception & Error_In) {
        } // for(b = 0; b < Num_Bodies; b++) {
      } // if(Print_Net_External_Forces == true) {
    } // if(t%TimeSteps_Between_Prints == 0) {

    #pragma omp single nowait
    {
      Print_time += Time_Since(time2);
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Apply Boundary conditions

    #pragma omp single nowait
    {
      time2 = Get_Time();
    } // #pragma omp single nowait

    #pragma omp single
    {
    for(b = 0; b < Num_Bodies; b++) {
      if(Bodies[b].Get_Is_Box() == true) { Apply_Box_BCs(Bodies[b], Box_Parameters[b]); }
    } // for(b = 0; b < Num_Bodies; b++)
    } // #pragma omp single

    #pragma omp single nowait
    {
      update_BC_time += Time_Since(time2);
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update Stress tensor (P)

    #pragma omp single nowait
    {
      time2 = Get_Time();
    } // #pragma omp single nowait

    for(b = 0; b < Num_Bodies; b++) {
      // Note: We don't update P for Bodys that are boundaries
      if(Bodies[b].Get_Boundary() == true) { continue; }

      else {
        /* Update each Particles's P tensor.
        We only update P when the bth Body's counter is zero.

        Note: the Update_P method has an orphaned for loop (and takes care of
        removing damaged particles in parallel, damaged particles are not
        removed until every particle's P tensor has been updated. This makes
        the code parallelizable and determinstic) */
        if(Time_Step_Index[b] == 0) {
          unsigned time_update = dt*Bodies[b].Get_Time_Steps_Between_Updates();
          Bodies[b].Update_P(time_update);
        } // if(Time_Step_Index[b] == 0) {
      } // else
    } // for(b = 0; b < Num_Bodies; b++) {

    #pragma omp single nowait
    {
      update_P_time += Time_Since(time2);
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
      time2 = Get_Time();
    } // #pragma omp single nowait

    /* First, we need to set each particle's contact force to zero. It should be
    noted that we only do this for a particular Body if that body
    is updating it's position this cycle. Otherwise, since the force won't be
    used for anything, there's no reason to waste CPU cycles setting that
    bodies's particle's contact forces to zero. */
    for(b = 0; b < Num_Bodies; b++) {
      unsigned Num_Particles = Bodies[b].Get_Num_Particles();

      #pragma omp for
      for(p = 0; p < Num_Particles; p++) {
        (Bodies[b])[p].Force_Contact = {0,0,0};
        (Bodies[b])[p].Force_Friction = {0,0,0};
      } // for(p = 0; p < (Bodies[b]).Get_Num_Particles(); p++) {
    } // for(b = 0; b < Num_Bodies; b++) {

    /* Now we can apply the contact algorythm. Note that this must be applied
    every time step no matter what (so that bodies that update each step can
    are proprly updated/have the right forces applied each timestpe) */
    for(unsigned b1 = 0; b1 < Num_Bodies - 1; b1++) {
      for(unsigned b2 = b1 + 1; b2 < Num_Bodies; b2++) {
        Body::Contact(Bodies[b2], Bodies[b1]);
      } // for(unsigned b2 = b1 + 1; b2 < Num_Bodies; b2++) {
    } // for(unsinged b1 = 0; b1 < Num_Bodies - 1; b1++) {

    #pragma omp single nowait
    {
      contact_time += Time_Since(time2);
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update Position (x)

    #pragma omp single nowait
    {
      time2 = Get_Time();
    } // #pragma omp single nowait

    for(b = 0; b < Num_Bodies; b++) {
      // Note: we don't update P for Bodies that are boundaries
      if(Bodies[b].Get_Boundary() == true) { continue; }

      else {
        /* We only want to update x (the traditional way) if we're on a timestep
        where the bth Body gets updated. Suppose that the bth body
        only updates once every k steps (meaning that Stpes_Between_Update[b] = k)
        on the 0th step, the bth Bodies's counter is zero. After
        each step it increments. On the kth step, its counter reaches k and the
        counter gets truncaed back to zero. Therefore, every k steps the bth
        Body's counter will be zero. Thus, we use a 0 counter
        as an indicator that we should update this Body. */
        if(Time_Step_Index[b] == 0) {
          /* First, update the 'F_Index' for the current Body. This
          controls which member of each particle's 'F' array is the 'newest'. */
          #pragma omp single
            Bodies[b].Increment_F_Index();

          // Now update the position of each particle in this body.
          Bodies[b].Update_x(dt);
        } // if(Time_Step_Index[b] == 0) {
        else {
          /* If we're not on an update step, then we'll let this body continue
          accelerating at whatever acceleration it attained after the last
          time step. */
          unsigned Num_Particles = (Bodies[b]).Get_Num_Particles();

          #pragma omp for
          for(p = 0; p < Num_Particles; p++) {
            if((Bodies[b])[p].Get_D() >= 1) { continue; }

            (Bodies[b])[p].x += dt*(Bodies[b])[p].V;       // x_p+1 = x_p + dt*v_(p+1/2)           : mm Vector
            (Bodies[b])[p].V += dt*(Bodies[b])[p].a;       // V_p+3/2 = V_p+1/2 + dt*a(t_p+1)      : mm/s Vector
          } // for(p = 0; p < (Bodies[b]).Get_Num_Particles(); p++) {
        } // else {
      } // else {
    } // for(b = 0; b < Num_Bodies; b++) {

    #pragma omp single nowait
    {
      update_x_time += Time_Since(time2);
    } // #pragma omp single nowait



    ////////////////////////////////////////////////////////////////////////////
    // Update each time step counter
    /* Here we increment each Body's counter. If a particular counter
    reaches its limit (the value of Bodies[b].Time_Steps_Between_Updates) then we set that
    counter to zero (reset the counter). */

    #pragma omp single
    {
    for(b = 0; b < Num_Bodies; b++) {
      Time_Step_Index[b]++;

      if(Time_Step_Index[b] == Bodies[b].Get_Time_Steps_Between_Updates()) {
        Time_Step_Index[b] = 0;
      } // if(Time_Step_Index[b] == Time_Setps_Between_Updates[b]) {
    } // for(b = 0; b < Num_Bodies; b++) {
    } // #pragma omp single
  } // for(t = 0; t < Num_Time_Steps; t++) {
  } // #pragma omp parallel
  printf(         "Done!\n");
  simulation_time = Time_Since(time1);

  // If saving is enabled, Dump particle data to file
  if(Save_Simulation == 1) {
    Data_Dump::Save_Simulation(Bodies, Num_Bodies);
  } // if(Save_Simulation == 1) {

  // Print timing data
  #if defined(_OPENMP)
    printf(       "\nIt took %lf s to perform %u Particle time steps \n",simulation_time, Num_Time_Steps);
    printf(       "%lfs to update BC's\n", update_BC_time);
    printf(       "%lfs to update P\n", update_P_time);
    printf(       "%lfs for Contact\n", contact_time);
    printf(       "%lfs to update x\n", update_x_time);
    printf(       "%lfs to print data to files\n", Print_time);

  #else
    unsigned long MS_Iter,                                           // These are used to store the number of miliseconds that
                  MS_BC,                                             // it took to execute each of the major operations in
                  MS_P,                                              // the code. These are only used the the code is executed
                  MS_Contact,                                        // sequentially
                  MS_x,
                  MS_Print;

    MS_Iter = (unsigned long)((double)time1 / (double)CLOCKS_PER_MS);
    MS_BC = (unsigned long)((double)update_BC_time / (double)CLOCKS_PER_MS);
    MS_P = (unsigned long)((double)update_P_time / (double)CLOCKS_PER_MS);
    MS_Contact = (unsigned long)((double)contact_time / (double)CLOCKS_PER_MS);
    MS_x = (unsigned long)((double)update_x_time / (double)CLOCKS_PER_MS);
    MS_Print = (unsigned long)((double)Print_time / (double)CLOCKS_PER_MS);

    printf(         "\nIt took %lu ms to perform %u Particle time steps \n",MS_Iter, Num_Time_Steps);
    printf(         "%lums to update BC's\n", MS_BC);
    printf(         "%lums to update P\n", MS_P);
    printf(         "%lums for Contact\n", MS_Contact);
    printf(         "%lums to update x\n", MS_x);
    printf(         "%lums to print data to files\n", MS_Print);
  #endif

  delete [] Bodies;
} // void Simulation(void) {



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
    Data_Dump::Load_Simulation(Bodies, Num_Bodies);

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
