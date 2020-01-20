#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "Load_Simulation.h"
#include "IO_Ops.h"
#include "Errors.h"
#include <string>
#include <stdio.h>

void IO::Load_Simulation(Body ** Bodies_Ptr, unsigned & Num_Bodies) {
  // First, open up the Particle_data file
  FILE * File = fopen("./IO/Saves/Simulation_Data.txt","r");
  if(File == nullptr) {
    throw Cant_Open_File("Can't Open File Exception: Thrown by IO::Load_Simulation\n"
                         "For some reason, /IO/Saves/Simulation_Data.txt could not be opened :(\n");
  } // if(File == nullptr) {

  fseek(File, 0, SEEK_SET);

  // Buffers
  unsigned Buf_Length = 100;
  char Buf[Buf_Length];                                 // Buffer to store text from file
  unsigned uBuf;
  std::string strBuf;

  // First, read how many Bodies we need to make.
  strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str(), " %u \n \n", &Num_Bodies);

  // Now use this information to make the Bodies array
  *Bodies_Ptr = new Body[Num_Bodies];

  // Now read in each Body's name and Fundamental properties
  unsigned Dimensions[3];
  unsigned Is_Box;
  for(unsigned i = 0; i < Num_Bodies; i++) {
    /* We want to get each particle's name. The issue here is that there's no
    way to directly scan to a string. Instead, we can scan to a char array, Buf
    in this case. This copies over the contents of the name as well as a null
    terminating character into Buf. We can then assign Buf to a string buffer.
    This works by copying over the characters from Buf into the string until
    the null terminating character is encountered. To ensure that this happens
    (which may not be the case if the string is too long), we manually assign
    the final element of Buf to \0 before doing this (Otherwise, we'd get a segmenetation
    fault if the name was longer than 100 characters). We then allocate a new
    body and assign its name to the string. */

    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str(), " %s \n", Buf);
    Buf[Buf_Length-1] = '\0';          // Incase sscanf filled Buf
    strBuf = Buf;
    #if defined(LOAD_MONITOR)
      printf("\nRead Body %u's name as                  %s\n", i, Buf);
    #endif


    // Set Bodys name.
    (*Bodies_Ptr)[i].Set_Name(strBuf);


    // Now determine if this Body is a Box
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &Is_Box);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Box as                %u\n", i, Is_Box);
    #endif


    if(Is_Box == true) {
      strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &Dimensions[0]);
      strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &Dimensions[1]);
      strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &Dimensions[2]);
      (*Bodies_Ptr)[i].Set_Box_Dimensions(Dimensions[0], Dimensions[1], Dimensions[2]);

      #if defined(LOAD_MONITOR)
        printf("Read Body %u's dimensions as            <%u %u %u>\n", i, Dimensions[0], Dimensions[1], Dimensions[2]);
      #endif
    } // if(uBuf == true) {


    // Now read in the 'Is fixed in place' flag
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Set_Is_Fixed(uBuf);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Fixed as              %u\n", i, uBuf);
    #endif


    // Now read in the 'Is damageable' flag
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Set_Damageable(uBuf);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Damageable as         %u\n", i, uBuf);
    #endif


    // Now read in number of particles and use if this Body is not a Box
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    if(Is_Box == false) { (*Bodies_Ptr)[i].Set_Num_Particles(uBuf); }
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Num Particles as         %u\n", i, uBuf);
    #endif


    // Now read in time steps between updates
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    if(Is_Box == false) { (*Bodies_Ptr)[i].Set_Time_Steps_Between_Updates(uBuf); }
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's time steps per update as %u\n", i, uBuf);
    #endif



    // Finally read in File number information
    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Net_External_Force = uBuf;

    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Particle_Forces = uBuf;

    strBuf = read_line_after_char(File, ':'); sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Particle_Positions = uBuf;
  } // for(unsigned i = 0; i < Num_Bodiess; i++) {

  // pass the newly created body's to the 'load body' function
  for(unsigned i = 0; i < Num_Bodies; i++) { IO::Load_Body((*Bodies_Ptr)[i]); }

  fclose(File);
} // void IO::Load_Simulation(Body ** Bodies_Ptr, unsigned & Num_Bodies) {



void IO::Load_Body(Body & Body_In) {
  /* This function is designed to read in particle and use it to create a Body */

  // First, open up the desired file.
  unsigned Num_Particles = 0;

  // First, get a path to the file
  std::string File_Path = "./IO/Saves/";
  File_Path += Body_In.Get_Name();
  File_Path += ".txt";

  FILE * File = fopen(File_Path.c_str(), "r");
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Can't Open File Exception: Thrown by IO::Load_Simulation\n"
            "For some reason, /IO/Saves/%s.txt could not be opened :(\n",
            Body_In.Get_Name().c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {

  // Buffers to hold variables that we read in (we need to do this b/c the
  // Body's varialbes are hidden/seed to be set with setters)
  unsigned Buf_Length = 100;
  char Buf[Buf_Length];                          // Buffer to store text from file
  unsigned uBuf;
  double lfBuf;
  Materials::Material Mat;


  // We already have the Body's name, Box/Fixed flags, and
  // dimensions (if Box). We can, therefore, skip over these lines.
  fgets(Buf, 99, File);                          // Skip 'name' line.
  fgets(Buf, 99, File);                          // Skip blank line
  fgets(Buf, 99, File);                          // Skip 'Is a Box' line

  if(Body_In.Get_Is_Box() == true) {                // if a Box, skip the 3 dimension lines.
    fgets(Buf, 99, File);
    fgets(Buf, 99, File);
    fgets(Buf, 99, File);
  } // if(Body_In.Get_Is_Box() == true) {

  fgets(Buf, 99, File);                          // Skip 'is fixed in place' line
  fgets(Buf, 99, File);                          // Skip 'Is Damageable' line
  fgets(Buf, 99, File);                          // Skip blank line
  fgets(Buf, 99, File);                          // Skip 'Kerenel-Parameters' line

  //////////////////////////////////////////////////////////////////////////////
  // Read in Kernel parameters.
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Body_In.Set_Inter_Particle_Spacing(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %u\n", &uBuf);      Body_In.Set_Support_Radius(uBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);         // Skip 'Support radius (mm)' line
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n\n", &lfBuf);       // Skip 'Shape Function Amplitude' line

  fgets(Buf, 99, File);                          // Skip 'Material-Parameters' line.

  //////////////////////////////////////////////////////////////////////////////
  // Read in Material parameters.

  // Read in material name
  fread(Buf, 1, 30, File); fscanf(File, " %s\n", Buf);
  Buf[Buf_Length-1] = '\0';
  Mat.Name = Buf;

  // Read in other material properties
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Mat.Lame = lfBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Mat.mu0 = lfBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Body_In.Set_mu(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %u\n", &uBuf);      Body_In.Set_F_Index(uBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Mat.E = lfBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Mat.density = lfBuf;

  // Our material is now fully characterized, we can set Particle's material
  Body_In.Set_Material(Mat);

  // Read other material properties.
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Body_In.Set_alpha(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n\n", &lfBuf);  Body_In.Set_Tau(lfBuf);

  //////////////////////////////////////////////////////////////////////////////
  // Read in Particle properties.

  // Now read in number of particles
  fgets(Buf, 99, File);                          // Skip 'Particles' line.
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Num_Particles);

  // Now read in particles.
  for(unsigned i = 0; i < Num_Particles; i++) { IO::Load_Particle(Body_In[i], File); }

  //////////////////////////////////////////////////////////////////////////////
  // Now recreate each Particle's neighbor arrays.
  unsigned Neighbor_ID;
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,                                // Shape Tensor (zero initialized)      : unitless Tensor
           0,0,0,
           0,0,0};
  const double Shape_Function_Amp = Body_In.Get_Shape_Function_Amplitude();
  const double h = Body_In.Get_h();

  for(unsigned i = 0; i < Num_Particles; i++) {
    // Check that the current particle has Neighbors

    if(Body_In[i].Num_Neighbors != 0) {
      // If so, then set up this particle's neighbor arrays.
      A = Tensor(0,0,0,
                 0,0,0,
                 0,0,0);

      for(unsigned j = 0; j < Body_In[i].Num_Neighbors; j++) {
        Neighbor_ID = Body_In[i].Neighbor_IDs[j];

        // Calculate displacement vectors
        Body_In[i].R[j] = Body_In[Neighbor_ID].X - Body_In[i].X;
        Body_In[i].Mag_R[j] = Body_In[i].R[j].Magnitude();

        // Calculate shape function, shape function gradient for jth neighbor
        Body_In[i].W[j] = Shape_Function_Amp*(h - Body_In[i].Mag_R[j])
                            *(h - Body_In[i].Mag_R[j])
                            *(h - Body_In[i].Mag_R[j]);

        Body_In[i].Grad_W[j] = (-3*Shape_Function_Amp
                                 *((h - Body_In[i].Mag_R[j])*(h - Body_In[i].Mag_R[j]))/ Body_In[i].Mag_R[j])
                                 *Body_In[i].R[j];

        // Add in the Current Neighbor's contribution to the Shape tensor
        V_j = Body_In[Neighbor_ID].Vol;
        A += Dyadic_Product((V_j*Body_In[i].Grad_W[j]), Body_In[i].R[j]);
      } // for(unsigned j = 0; j < Body_In[i].Num_Neighbors; i++) {

      // Now we can calculate A^(-1) from A.
      Body_In[i].A_Inv = A^(-1);

      // Now that neighbors have been set, we set 'Neighbors_Are_Set' to true
      Body_In[i].Neighbors_Are_Set = true;
    } // if(Body_In[i].Num_Neighbors != 0) {
  } // for(unsigned i = 0; i < Body_In; i++) {

  // All done, close the file.
  fclose(File);
} // void IO::Load_Body(Body & Body_In) {



void IO::Load_Particle(Particle & P_In, FILE * File) {
  /* This function reads in the particle data for a specific particle. This
  data is transfered from the File to P_In. */

  char Buf[100];                                 // Buffer to store text from file (discarded)

  // First, read in the particle's ID and dimensions
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &P_In.ID);

  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Mass);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Vol);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Radius);

  // Now read in particle dynamic properties
  fread(Buf, 1, 30, File); fscanf(File, " <%lf %lf %lf>\n", &P_In.X(0), &P_In.X(1), &P_In.X(2));
  fread(Buf, 1, 30, File); fscanf(File, " <%lf %lf %lf>\n", &P_In.x(0), &P_In.x(1), &P_In.x(2));
  fread(Buf, 1, 30, File); fscanf(File, " <%lf %lf %lf>\n", &P_In.V(0), &P_In.V(1), &P_In.V(2));
  fread(Buf, 1, 30, File); fscanf(File, " |%lf %lf %lf|\n", &P_In.F[0](0,0), &P_In.F[0](0,1), &P_In.F[0](0,2));
                           fscanf(File, " |%lf %lf %lf|\n", &P_In.F[0](1,0), &P_In.F[0](1,1), &P_In.F[0](1,2));
                           fscanf(File, " |%lf %lf %lf|\n", &P_In.F[0](2,0), &P_In.F[0](2,1), &P_In.F[0](2,2));
  fread(Buf, 1, 30, File); fscanf(File, " |%lf %lf %lf|\n", &P_In.F[1](0,0), &P_In.F[1](0,1), &P_In.F[1](0,2));
                           fscanf(File, " |%lf %lf %lf|\n", &P_In.F[1](1,0), &P_In.F[1](1,1), &P_In.F[1](1,2));
                           fscanf(File, " |%lf %lf %lf|\n", &P_In.F[1](2,0), &P_In.F[1](2,1), &P_In.F[1](2,2));

  // Damage parameters
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_H);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_M);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_Critical);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.D);

  // Neighbor paramaters
  P_In.Neighbors_Are_Set = false;
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &P_In.Num_Neighbors);

  // Now allocate memory for P_In's neighbor arrays
  P_In.Neighbor_IDs = new unsigned[P_In.Num_Neighbors];
  P_In.R = new Vector[P_In.Num_Neighbors];                                     //        : mm Vector
  P_In.Mag_R = new double[P_In.Num_Neighbors];                                 //        : mm
  P_In.W = new double[P_In.Num_Neighbors];                                     //        : unitless
  P_In.Grad_W = new Vector[P_In.Num_Neighbors];                                //        : 1/mm Vector

  // Now read in neighbor IDs. Before we can do that, however, we need to move
  // the file pointer ahead, past 'Neighbor IDs: '
  fread(Buf, 1, 30, File);
  for(unsigned i = 0; i < P_In.Num_Neighbors; i++) {
    fscanf(File, " %u", &P_In.Neighbor_IDs[i]);
  } // for(unsigned i = 0; i < P_In.Num_Neighbors; i++) {
} // void IO::Load_Particle(Particle & P_In, FILE * File) {
