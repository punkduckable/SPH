#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Tensor/Tensor.h"
#include "Vector/Vector.h"
#include "Load_Simulation.h"
#include "IO_Ops.h"
#include "Errors.h"
#include <string>
#include <fstream>
#include <stdio.h>

void IO::Load_Simulation(Body ** Bodies_Ptr, unsigned & Num_Bodies) {
  // First, open up the Particle_data file
  std::ifstream File;
  File.open("./IO/Saves/Simulation_Data.txt");
  if(File.is_open() == false) {
    throw Cant_Open_File("Can't Open File Exception: Thrown by IO::Load_Simulation\n"
                         "For some reason, /IO/Saves/Simulation_Data.txt could not be opened :(\n");
  } // if(File.is_open() == false) {

  // Buffers
  unsigned uBuf;
  std::string strBuf;

  // First, read how many Bodies we need to make.
  strBuf = read_line_after(File, "Number of Bodies:");
  sscanf(strBuf.c_str(), " %u \n \n", &Num_Bodies);

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

    strBuf = read_line_after(File, "name:");
    unsigned Buf_Length = 256;
    char Buf[Buf_Length];
    sscanf(strBuf.c_str(), " %s \n", Buf);
    Buf[Buf_Length-1] = '\0';          // In case sscanf filled Buf
    strBuf = Buf;
    #if defined(LOAD_MONITOR)
      printf("\nRead Body %u's name as                  %s\n", i, Buf);
    #endif


    // Set Bodys name.
    (*Bodies_Ptr)[i].Set_Name(strBuf);


    // Now determine if this Body is a Box
    strBuf = read_line_after(File, "Is a Box:");
    sscanf(strBuf.c_str()," %u \n", &Is_Box);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Box as                %u\n", i, Is_Box);
    #endif


    if(Is_Box == true) {
      strBuf = read_line_after(File, "X_SIDE_LENGTH:"); sscanf(strBuf.c_str()," %u \n", &Dimensions[0]);
      strBuf = read_line_after(File, "Y_SIDE_LENGTH:"); sscanf(strBuf.c_str()," %u \n", &Dimensions[1]);
      strBuf = read_line_after(File, "Z_SIDE_LENGTH:"); sscanf(strBuf.c_str()," %u \n", &Dimensions[2]);
      (*Bodies_Ptr)[i].Set_Box_Dimensions(Dimensions[0], Dimensions[1], Dimensions[2]);

      #if defined(LOAD_MONITOR)
        printf("Read Body %u's dimensions as            <%u %u %u>\n", i, Dimensions[0], Dimensions[1], Dimensions[2]);
      #endif
    } // if(uBuf == true) {


    // Now read in the 'Is fixed in place' flag
    strBuf = read_line_after(File, "Is fixed in place:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Set_Is_Fixed(uBuf);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Fixed as              %u\n", i, uBuf);
    #endif


    // Now read in the 'Is damageable' flag
    strBuf = read_line_after(File, "Is Damageable:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Set_Damageable(uBuf);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Is_Damageable as         %u\n", i, uBuf);
    #endif


    // Now read in number of particles and use if this Body is not a Box
    strBuf = read_line_after(File, "Number of particles:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    if(Is_Box == false) { (*Bodies_Ptr)[i].Set_Num_Particles(uBuf); }
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's Num Particles as         %u\n", i, uBuf);
    #endif


    // Now read in time steps between updates
    strBuf = read_line_after(File, "Time steps per update:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Set_Time_Steps_Between_Updates(uBuf);
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's time steps per update as %u\n", i, uBuf);
    #endif



    // Finally read in File number information
    strBuf = read_line_after(File, "# times printed net external forces:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Net_External_Force = uBuf;
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's # times printed net external forces as %u\n", i, uBuf);
    #endif

    strBuf = read_line_after(File, "# times printed particle forces:");
    sscanf(strBuf.c_str()," %u \n \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Particle_Forces = uBuf;
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's # times printed particle forces as %u\n", i, uBuf);
    #endif

    strBuf = read_line_after(File, "# times printed particle positions:");
    sscanf(strBuf.c_str()," %u \n", &uBuf);
    (*Bodies_Ptr)[i].Times_Printed_Particle_Positions = uBuf;
    #if defined(LOAD_MONITOR)
      printf("Read Body %u's # times printed particle positions as %u\n", i, uBuf);
    #endif
  } // for(unsigned i = 0; i < Num_Bodiess; i++) {

  // pass the newly created body's to the 'load body' function
  for(unsigned i = 0; i < Num_Bodies; i++) { IO::Load_Body((*Bodies_Ptr)[i]); }

  File.close();
} // void IO::Load_Simulation(Body ** Bodies_Ptr, unsigned & Num_Bodies) {



void IO::Load_Body(Body & Body_In) {
  /* This function is designed to read in particle and use it to create a Body */

  // First, open up the desired file.
  std::string File_Path = "./IO/Saves/";
  File_Path += Body_In.Get_Name();
  File_Path += ".txt";

  std::ifstream File;
  File.open(File_Path);

  if(File.is_open() == false) {
    char Buf[500];
    sprintf(Buf,
            "Can't Open File Exception: Thrown by IO::Load_Simulation\n"
            "For some reason, %s could not be opened :(\n",
            File_Path.c_str());
    throw Cant_Open_File(Buf);
  } // if(File.is_open() == false) {


  // Buffers
  std::string strBuf;
  unsigned uBuf;
  double lfBuf;

  /* We already know the Body's name, Box/Fixed flaes, and dimensions (if
  the body is a box). These values were read in by IO::Load_Simulation.
  Therefore, we can skip over these lines. */
  IO::read_line_after(File, "Name: ");
  IO::read_line_after(File, "Is a Box: ");
  if(Body_In.Get_Is_Box() == true) {
    IO::read_line_after(File, "X_SIDE_LENGTH: ");
    IO::read_line_after(File, "Y_SIDE_LENGTH: ");
    IO::read_line_after(File, "Z_SIDE_LENGTH: ");
  }
  IO::read_line_after(File, "Is Fixed in place: ");
  IO::read_line_after(File, "Is Damageable: ");



  //////////////////////////////////////////////////////////////////////////////
  // Read in Kernel parameters.
  strBuf = IO::read_line_after(File, "Inter Particle Spacing:");
  sscanf(strBuf.c_str(), " %lf \n", &lfBuf);
  Body_In.Set_Inter_Particle_Spacing(lfBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's IPS as:                       %lf\n",Body_In.Get_Name().c_str(), lfBuf);
  #endif

  strBuf = IO::read_line_after(File, "Support Radius (IPS):");
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  Body_In.Set_Support_Radius(uBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's Support Radius as:            %u\n",Body_In.Get_Name().c_str(), uBuf);
  #endif

  // Skip these lines (they're redundant)
  IO::read_line_after(File, "Support Radius (mm) aka h:");
  IO::read_line_after(File, "Shape Function Amplitude:");



  //////////////////////////////////////////////////////////////////////////////
  // Read in Material parameters.

  unsigned Buf_Length = 256;
  char Buf[Buf_Length];
  Materials::Material Mat;

  // Read in material name
  strBuf = IO::read_line_after(File, "Material:");
  sscanf(strBuf.c_str()," %s \n", Buf);
  Buf[Buf_Length-1] = '\0';                      // This ensures that Buf is a null-terminated string
  Mat.Name = Buf;
  #if defined(LOAD_MONITOR)
    printf("Read %s's Material name as:             %s\n", Body_In.Get_Name().c_str(), Buf);
  #endif

  // Read in Lame parameter
  strBuf = IO::read_line_after(File, "Lame parameter:");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Mat.Lame = lfBuf;
  #if defined(LOAD_MONITOR)
    printf("Read %s's Lame parameter as:            %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Read in Shear modulus (mu0)
  strBuf = IO::read_line_after(File, "Shear modulus (mu0):");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Mat.mu0 = lfBuf;
  #if defined(LOAD_MONITOR)
    printf("Read %s's Shear modulus (mu0) as:       %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Read in Viscosity (mu)
  strBuf = IO::read_line_after(File, "Viscosity (mu):");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Body_In.Set_mu(lfBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's Viscosity (mu) as:            %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Read in F_Index
  strBuf = IO::read_line_after(File, "F_Index:");
  sscanf(strBuf.c_str()," %u \n", &uBuf);
  Body_In.Set_F_Index(uBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's F_Index as %u\n", Body_In.Get_Name().c_str(), uBuf);
  #endif

  // Read in Hourglass Stiffness (E)
  strBuf = IO::read_line_after(File, "Hourglass Stiffness (E):");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Mat.E = lfBuf;
  #if defined(LOAD_MONITOR)
    printf("Read %s's Hourglass Stuffness (E) as:   %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Read in Material density
  strBuf = IO::read_line_after(File, "Material density:");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Mat.density = lfBuf;
  #if defined(LOAD_MONITOR)
    printf("Read %s's density as:                   %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Our material is now fully characterized, we can set Particle's material
  Body_In.Set_Material(Mat);

  // Read in alpha (HG parameter)
  strBuf = IO::read_line_after(File, "alpha (HG parameter):");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Body_In.Set_alpha(lfBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's alpha (HG parameter) as:      %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  // Read in Tau (damage parameter)
  strBuf = IO::read_line_after(File, "Tau (damage parameter):");
  sscanf(strBuf.c_str()," %lf \n", &lfBuf);
  Body_In.Set_Tau(lfBuf);
  #if defined(LOAD_MONITOR)
    printf("Read %s's Tau (damage parameter) as:    %lf\n", Body_In.Get_Name().c_str(), lfBuf);
  #endif

  //////////////////////////////////////////////////////////////////////////////
  // Read in Particle properties.

  // read in number of particles
  unsigned Num_Particles;
  strBuf = IO::read_line_after(File, "Number of particles:");
  sscanf(strBuf.c_str(), " %u \n", &Num_Particles);
  #if defined(LOAD_MONITOR)
    printf("Read %s's Number of particles as:       %u\n", Body_In.Get_Name().c_str(), Num_Particles);
  #endif

  /* Now read in particles. For each particle, this sets the particles members
  that do not depend on the particles neighbors. The neighbor dependnet
  members are set up later. The reason for this is that the particles are
  read in one-by-one, meaning that neighbor information will not be available
  until every particle has been read in. */
  for(unsigned i = 0; i < Num_Particles; i++) { IO::Load_Particle(Body_In[i], File); }

  /* Now, for each particle in the Body, set the particles neighbor dependent
  members. */
  for(unsigned i = 0; i < Num_Particles; i++) { Body_In.Set_Neighbor_Dependent_Members(i); }

  // All done, close the file.
  File.close();
} // void IO::Load_Body(Body & Body_In) {



void IO::Load_Particle(Particle & P_In, std::ifstream & File) {
  /* This function reads in the particle data for a specific particle. This
  data is transfered from the File to P_In. */

  char Buf[100];                                 // Buffer to store text from file (discarded)
  std::string strBuf;                            // Another buffer to store text


  //////////////////////////////////////////////////////////////////////////////
  // Particle's ID and dimensions

  strBuf = read_line_after(File, "ID:");         sscanf(strBuf.c_str(), " %u \n", &P_In.ID);

  strBuf = read_line_after(File, "Mass:");       sscanf(strBuf.c_str(), " %lf \n", &P_In.Mass);
  strBuf = read_line_after(File, "Volume:");     sscanf(strBuf.c_str(), " %lf \n", &P_In.Volume);
  strBuf = read_line_after(File, "Radius:");     sscanf(strBuf.c_str(), " %lf \n", &P_In.Radius);


  //////////////////////////////////////////////////////////////////////////////
  // Particle dynamic properties

  strBuf = read_line_after(File, "X:");         sscanf(strBuf.c_str(), " <%lf %lf %lf> \n", &P_In.X(0), &P_In.X(1), &P_In.X(2));
  strBuf = read_line_after(File, "x:");         sscanf(strBuf.c_str(), " <%lf %lf %lf> \n", &P_In.x(0), &P_In.x(1), &P_In.x(2));
  strBuf = read_line_after(File, "V:");         sscanf(strBuf.c_str(), " <%lf %lf %lf> \n", &P_In.V(0), &P_In.V(1), &P_In.V(2));

  strBuf = read_line_after(File, "F[0]:");      sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[0](0,0), &P_In.F[0](0,1), &P_In.F[0](0,2));
  strBuf = read_line_after(File, " ");          sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[0](1,0), &P_In.F[0](1,1), &P_In.F[0](1,2));
  strBuf = read_line_after(File, " ");          sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[0](2,0), &P_In.F[0](2,1), &P_In.F[0](2,2));

  strBuf = read_line_after(File, "F[1]:");      sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[1](0,0), &P_In.F[1](0,1), &P_In.F[1](0,2));
  strBuf = read_line_after(File, " ");          sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[1](1,0), &P_In.F[1](1,1), &P_In.F[1](1,2));
  strBuf = read_line_after(File, " ");          sscanf(strBuf.c_str(), " |%lf %lf %lf| \n", &P_In.F[1](2,0), &P_In.F[1](2,1), &P_In.F[1](2,2));


  //////////////////////////////////////////////////////////////////////////////
  // Damage parameters

  strBuf = read_line_after(File, "Stretch_H:");  sscanf(strBuf.c_str(), " %lf \n", &P_In.Stretch_H);
  strBuf = read_line_after(File, "Stretch_M:");  sscanf(strBuf.c_str(), " %lf \n", &P_In.Stretch_M);
  strBuf = read_line_after(File, "Stretch_Critical:");     sscanf(strBuf.c_str(), " %lf \n", &P_In.Stretch_Critical);
  strBuf = read_line_after(File, "D:");          sscanf(strBuf.c_str(), " %lf \n", &P_In.D);


  //////////////////////////////////////////////////////////////////////////////
  // BC Information

  // We need more buffers to read in BCs
  char Buffers[3][10];

  // Read in "Has_BC", write it to P_In
  strBuf = read_line_after(File, "Has Boundary Conditions:");
  sscanf(strBuf.c_str(), " < %s %s %s > \n", Buffers[0], Buffers[1], Buffers[2]);
  for(unsigned i = 0; i < 3; i++) {
    if(String_Ops::Contains(Buffers[i], "true")) { P_In.Has_BC[i] = true; }
    else { P_In.Has_BC[i] = false; }
  } // for(unsigned i = 0; i < 3; i++) {

  // Read in "BC", write to P_In
  strBuf = read_line_after(File, "Boundary Conditions:");
  sscanf(strBuf.c_str(), " < %s %s %s > \n", Buffers[0], Buffers[1], Buffers[2]);
  for(unsigned i = 0; i < 3; i++) {
    if(P_In.Has_BC[i] == true ) { sscanf(Buffers[i]," %lf ", &P_In.BC[i]); }
  } // for(unsigned i = 0; i < 3; i++) {


  //////////////////////////////////////////////////////////////////////////////
  // Neighbor paramaters

  // Read in number of neighbors
  P_In.Neighbors_Are_Set = false;
  unsigned Num_Neighbors;
  strBuf = read_line_after(File, "Number of neighbors:");  sscanf(strBuf.c_str(), " %u \n", &Num_Neighbors);

  // Now read in neighbor IDs into a list. Before we can do that, however,
  // we need to move the file pointer ahead, past 'Neighbor IDs: '
  File.get(Buf, 30);
  P_In.Num_Neighbors = Num_Neighbors;

  // Set Particle's Neighbor IDs
  unsigned * Neighbor_IDs = new unsigned[Num_Neighbors];
  for(unsigned j = 0; j < Num_Neighbors; j++) { File >> Neighbor_IDs[j]; }
  P_In.Neighbor_IDs = Neighbor_IDs;
} // void IO::Load_Particle(Particle & P_In, std::ifstream & File) {
