#if !defined(DATA_DUMP_SOURCE)
#define DATA_DUMP_SOURCE

#include "Data_Dump.h"

void Data_Dump::Save_Simulation(const Body * Arrays, const unsigned Num_Arrays) {
  /* This Function is used to save a simulation. This function prints all the
  information that is needed to re-create the state of a simulation. This is
  done by first printing information on the number of particle arrays, as well
  as their names. Once this is done, each particle array gets a new file that
  contains all the information needed to re-create that Body object.
  The intent of this is to allow the user 'save' a simulation, allowing the user
  to - at a future point - pick up where they left off.*/

  unsigned i,j;

  // First, we need to know how many times each Body has been printed
  unsigned * VTK_File_Numbers = new unsigned[Num_Arrays];
  unsigned * Force_File_Numbers = new unsigned[Num_Arrays];

  for(i = 0;  i < Num_Arrays; i++) {
    // First, get the ith particle's array VTK file number
    for(j = 0; j < VTK_File::Name_List.Node_Count(); j++)
      if(VTK_File::Name_List[j].compare(Arrays[i].Get_Name()) == 0) {
        VTK_File_Numbers[i] = VTK_File::File_Number_List[j];
        break;
      } // if(VTK_File::Name_List[j].compare(Arrays[i].Get_Name()) == 0) {

    /* Now we need to check if the ith particle array is not in the list. This
    would happen if j == VTK_File::Name_List.Node_Count (if the for loop
    never broke). If this is the case, then we should set VTK_File_Number[i]
    to zero. */
    if(j == VTK_File::Name_List.Node_Count())
        VTK_File_Numbers[i] = 0;


    // Now, get the ith particle's array Force file number
    for(j = 0; j < Particle_Debugger::Name_List.Node_Count(); j++)
      if(Particle_Debugger::Name_List[j].compare(Arrays[i].Get_Name()) == 0) {
        Force_File_Numbers[i] = Particle_Debugger::File_Number_List[j];
        break;
      } // if(Particle_Debugger::Name_List[j].compare(Arrays[i].Get_Name()) == 0) {

    /* Now we need to check if the ith particle array is not in the list. This
    would happen if j == Particle_Debugger::Name_List.Node_Count (if the for loop
    never broke). If this is the case, then we should set Force_File_Number[i]
    to zero. */
    if(j == Particle_Debugger::Name_List.Node_Count())
        Force_File_Numbers[i] = 0;
  } // for(i = 0;  i < Num_Arrays; i++) {

  // Open a new file. This file will store the number of Bodys as
  // well as their names.
  FILE * File = fopen("../Files/Saves/Body_Data.txt","w");

  // Print number of Arrays to this file
  fprintf(File,   "Number of Particle Arrays:    %u\n\n", Num_Arrays);

  // Now print each Body's name and other essential information
  for(i = 0; i < Num_Arrays; i++) {
    fprintf(File, "Body %3u name:      %s\n", i, Arrays[i].Get_Name().c_str());
    fprintf(File, "     Is a Cuboid:             %u\n",    Arrays[i].Get_Cuboid());

    if(Arrays[i].Get_Cuboid() == true) {
      fprintf(File, "          X_SIDE_LENGTH:      %u\n",  Arrays[i].Get_X_SIDE_LENGTH());
      fprintf(File, "          Y_SIDE_LENGTH:      %u\n",  Arrays[i].Get_Y_SIDE_LENGTH());
      fprintf(File, "          Z_SIDE_LENGTH:      %u\n",  Arrays[i].Get_Z_SIDE_LENGTH());
    } // if(Arrays[i].Get_Cuboid() == true) {

    fprintf(File, "     Is a Boundary:           %u\n",    Arrays[i].Get_Boundary());
    fprintf(File, "     Is Damageable:           %u\n",    Arrays[i].Get_Damagable());
    fprintf(File, "     Number of particles:     %u\n",    Arrays[i].Get_Num_Particles());

    fprintf(File, "     VTK_File number:         %u\n",    VTK_File_Numbers[i]);
    fprintf(File, "     Force_File number:       %u\n\n",  Force_File_Numbers[i]);
  } //   for(i = 0; i < Num_Arrays; i++) {

  // Now print each Particle array to its own file
  for(i = 0; i < Num_Arrays; i++)
    Save_Body(Arrays[i]);

  fclose(File);
} // void Data_Dump::Save_Simulation(const Body * Arrays, const unsigned Num_Arrays) {



void Data_Dump::Save_Body(const Body & Body_In) {
  /* This function prints all data needed to reconstruct a particular
  Body. This information is printed to a file named after the
  Body that it's printing from. This prints the basic material
  properties of the particles array, as well as all the information that is
  needed to recreate all the particles in that Body. The intent of
  this is to allow the user 'save' a Body object. */

  // Get name + number of particles of passed Body object.
  const std::string Name = Body_In.Get_Name();
  const unsigned Num_Particles = Body_In.Get_Num_Particles();

  // Now make a file path for this Particle's file
  std::string File_Path = "../Files/Saves/";
  File_Path += Name;
  File_Path += ".txt";

  /* Now let's make a file for this Body. Note, if there is already
  a save file for this Body it will be overwritten by this new file */
  FILE * File = fopen(File_Path.c_str(), "w");

  // Let's begin by printing the Body paramaters
  fprintf(File,   "Name:                         %s\n\n",  Name.c_str());

  fprintf(File,   "Is a cuboid:                  %u\n",    Body_In.Get_Cuboid());
  if(Body_In.Get_Cuboid() == true) {
    fprintf(File, "     X_SIDE_LENGTH:           %u\n",    Body_In.Get_X_SIDE_LENGTH());
    fprintf(File, "     Y_SIDE_LENGTH:           %u\n",    Body_In.Get_Y_SIDE_LENGTH());
    fprintf(File, "     Z_SIDE_LENGTH:           %u\n",    Body_In.Get_Z_SIDE_LENGTH());
  } //   if(Body_In.Get_Cuboid() == true) {
  fprintf(File,   "Is a Boundary:                %u\n",    Body_In.Get_Boundary());
  fprintf(File,   "Is Damageable:                %u\n\n",  Body_In.Get_Damagable());

  fprintf(File,   "       -- Kernel Parameters --\n");
  fprintf(File,   "Inter Particle Spacing:       %5lf\n",  Body_In.Get_Inter_Particle_Spacing());
  fprintf(File,   "Support Radius (IPS):         %u\n",    Body_In.Get_Support_Radius());
  fprintf(File,   "Support Radius (mm) aka h:    %5lf\n",  Body_In.Get_h());
  fprintf(File,   "Shape Function Amplitude:     %5lf\n\n",Body_In.Get_Shape_Function_Amplitude());

  fprintf(File,   "       -- Material Parameters --\n");
  fprintf(File,   "Material:                     %s\n",    Body_In.Get_Material().Name.c_str());
  fprintf(File,   "Lame parameter:               %5lf\n",  Body_In.Get_Lame());
  fprintf(File,   "Shear modulus (mu0):          %5lf\n",  Body_In.Get_mu0());
  fprintf(File,   "Viscosity (mu):               %5lf\n",  Body_In.Get_mu());
  fprintf(File,   "F_Index:                      %u\n",    Body_In.Get_F_Index());
  fprintf(File,   "Hourglass Stiffness (E):      %5lf\n",  Body_In.Get_E());
  fprintf(File,   "Material density:             %5lf\n",  Body_In.Get_density());
  fprintf(File,   "alpha (HG parameter):         %5lf\n",  Body_In.Get_alpha());
  fprintf(File,   "Tau (damage parameter):       %5lf\n\n",Body_In.Get_Tau());

  // Now let's print the number of particles
  fprintf(File,   "       -- Particles --\n");
  fprintf(File,   "Number of particles:          %u\n\n",    Body_In.Get_Num_Particles());

  // Finally, let's print the cuboid paramaters (should be removed if not using
  // a cuboid)
  //fprintf(File,   "X Side Length:                %u\n",    Simulation::X_SIDE_LENGTH);
  //fprintf(File,   "Y Side Length:                %u\n",    Simulation::Y_SIDE_LENGTH);
  //fprintf(File,   "Z Side Length:                %u\n\n",  Simulation::Z_SIDE_LENGTH);

  // Now let's print all particle data to the file
  for(unsigned i = 0; i < Num_Particles; i++)
    Save_Particle(Body_In[i], File, Body_In.Get_Cuboid());

  // We've now written the 'Particle_Data' file, we can close it.
  fclose(File);
} // void Data_Dump::Save_Body(const Body & Body_In) {



void Data_Dump::Save_Particle(const Particle & P_In, FILE * File, const bool Is_Cuboid) {
  /* This function prints all the information that is needed to re-create the
  input particle. Notably, this means that we do NOT need to print the first
  Piola Kirchoff stress tensor (P), the deformation gradient (F), any of the
  forces, or most of the neighbor arrays (the ID's are needed, the rest is not).
  P, F, and the forces are not needed because these are all calculated from
  stratch each iteration. Likewise, the neighbor array parameters can be
  recalculated if we know the neighbor IDs.

  Not storing this information in the File makes the file take up less/easier
  to read. This function assumes that the File has already been setup (with
  static particle class paramaters). */

  unsigned i;                                // index variable

  const Vector X = P_In.Get_X();
  const Vector x = P_In.Get_x();
  const Vector V = P_In.Get_V();
  const Tensor F_0 = P_In.Get_F(0);
  const Tensor F_1 = P_In.Get_F(1);

  // Print particle ID, dimensions
  fprintf(File,   "ID:                           %u\n",    P_In.Get_ID());
  fprintf(File,   "Mass:                         %5e\n",   P_In.Get_Mass());
  fprintf(File,   "Volume:                       %5e\n",   P_In.Get_Vol());
  fprintf(File,   "Radius:                       %5lf\n",  P_In.Get_Radius());

  // Print Particle dynamic properties
  fprintf(File,   "X:                            <%6.3lf %6.3lf %6.3lf>\n", X(0), X(1), X(2));
  fprintf(File,   "x:                            <%6.3lf %6.3lf %6.3lf>\n", x(0), x(1), x(2));
  fprintf(File,   "V:                            <%6.3lf %6.3lf %6.3lf>\n", V(0), V(1), V(2));
  fprintf(File,   "F[0]:                         |%6.3lf %6.3lf %6.3lf|\n", F_0(0,0), F_0(0,1), F_0(0,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_0(1,0), F_0(1,1), F_0(1,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_0(2,0), F_0(2,1), F_0(2,2));
  fprintf(File,   "F[1]:                         |%6.3lf %6.3lf %6.3lf|\n", F_1(0,0), F_1(0,1), F_1(0,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_1(1,0), F_1(1,1), F_1(1,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_1(2,0), F_1(2,1), F_1(2,2));


  // Damage paramaters
  fprintf(File,   "Stretch_H:                    %5lf\n",  P_In.Get_Stretch_H());
  fprintf(File,   "Stretch_M:                    %5lf\n",  P_In.Get_Stretch_M());
  fprintf(File,   "Stretch_Critical:             %5lf\n",  P_In.Get_Stretch_Critical());
  fprintf(File,   "D:                            %5lf\n",  P_In.Get_D());

  // Now, let's figure out how many neighbors this particle has.
  unsigned Num_Neighbors = P_In.Get_Num_Neighbors();

  // Neighbor paramters
  fprintf(File,   "Number of neighbors:          %u\n", P_In.Get_Num_Neighbors());

  // Print neighbor IDs
  fprintf(File,   "Neighbor IDs                  ");
  for(i = 0; i < Num_Neighbors; i++)
    fprintf(File,"%d ",P_In.Get_Neighbor_IDs(i));

  fprintf(File,"\n\n");
} // void Data_Dump::Save_Particle(const Particle & P_In, FILE * File, const bool Is_Cuboid) {



int Data_Dump::Load_Simulation(Body ** Array_Ptr, unsigned & Num_Arrays) {
  // First, open up the Particle_data file
  FILE * File = fopen("../Files/Saves/Body_Data.txt","r");
  fseek(File, 0, SEEK_SET);

  if(File == NULL) {
    printf("Couldn't find Body_Data.txt   :(\n");
    return 1;
  } // if(File == NULL) {

  // First, read how many Particle arrays we need to make.
  unsigned Buf_Length = 100;
  char Buf[Buf_Length];                                 // Buffer to store text from file
  unsigned uBuf;
  std::string strBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %u\n\n", &Num_Arrays);

  // Now use this information to genrate our arrays
  *Array_Ptr = new Body[Num_Arrays];

  Vector Dimensions;
  unsigned Is_Cuboid;

  // Now read in each array's name + Fundamental properties
  for(unsigned i = 0; i < Num_Arrays; i++) {
    /* We want to get each particle's name. The issue here is that there's no
    way to directly scan to a string. Instead, we can scan to a char array, Buf
    in this case. This copies over the contents of the name as well as a null
    terminating character into Buf. We can then assign Buf to a string buffer.
    This works by copying over the characters from Buf into the string until
    the null terminating character is encountered. To ensure that this happens
    (which may not be the case if the string is too long), we manually assign
    the final element of Buf to \0 before doing this (Otherwise, we'd get a segmenetation
    fault if the name was longer than 100 characters). We then allocate a new
    particle array and assign its name to the string. */

    fread(Buf, 1, 30, File); fscanf(File, " %s\n", Buf);
    Buf[Buf_Length-1] = '\0';
    strBuf = Buf;

    // Set Bodys name.
    (*Array_Ptr)[i].Set_Name(strBuf);

    // Now determine if this Body is a Cuboid
    fread(Buf, 1, 25, File); fscanf(File," %u\n", &Is_Cuboid);

    if(Is_Cuboid == true) {
      fread(Buf, 1, 20, File); fscanf(File," %u\n", &uBuf);   Dimensions(0) = uBuf;
      fread(Buf, 1, 20, File); fscanf(File," %u\n", &uBuf);   Dimensions(1) = uBuf;
      fread(Buf, 1, 20, File); fscanf(File," %u\n", &uBuf);   Dimensions(2) = uBuf;
      (*Array_Ptr)[i].Set_Cuboid_Dimensions(Dimensions);
    } // if(uBuf == true) {

    // Now read in the 'Is a boundary' flag
    fread(Buf, 1, 25, File); fscanf(File," %u\n", &uBuf);
    (*Array_Ptr)[i].Set_Boundary(uBuf);

    // Now read in the 'Is damageabl' flag
    fread(Buf, 1, 25, File); fscanf(File," %u\n", &uBuf);
    (*Array_Ptr)[i].Set_Damageable(uBuf);

    // Now read in number of particles and use if this Body is not a cuboid
    fread(Buf, 1, 25, File); fscanf(File," %u\n", &uBuf);
    if(Is_Cuboid == false) {
      (*Array_Ptr)[i].Set_Num_Particles(uBuf);
    } // if(Is_Cuboid == false) {

    // Finally read in File number information
    fread(Buf, 1, 25, File); fscanf(File," %u\n", &uBuf);
    VTK_File::Name_List.Add_Back((*Array_Ptr)[i].Get_Name());
    VTK_File::File_Number_List.Add_Back(uBuf);

    fread(Buf, 1, 25, File); fscanf(File," %u\n\n", &uBuf);
    Particle_Debugger::Name_List.Add_Back((*Array_Ptr)[i].Get_Name());
    Particle_Debugger::File_Number_List.Add_Back(uBuf);
  } // for(unsigned i = 0; i < Num_Arrays; i++) {

  // pass the newly created particle arrays to the 'load particle array' function
  for(unsigned i = 0; i < Num_Arrays; i++)
    Load_Body((*Array_Ptr)[i]);

  fclose(File);
  return 0;
} // int Data_Dump::Load_Saved_Simulation(Body ** Array_Ptr, unsigned & Num_Arrays) {



int Data_Dump::Load_Body(Body & Body_In) {
  /* This function is designed to read in particle and use it to create a
  particles array. */

  // First, open up the desired file.
  unsigned Num_Particles = 0;

  // First, get a path to the file
  std::string File_Path = "../Files/Saves/";
  File_Path += Body_In.Get_Name();
  File_Path += ".txt";

  FILE * File = fopen(File_Path.c_str(), "r");

  if(File == NULL) {
    printf("Could not open Particle data file\n");
    return 1;
  } // if(File == NULL) {

  // Buffers to hold variables that we read in (we need to do this b/c the
  // Body's varialbes are hidden/seed to be set with setters)
  unsigned Buf_Length = 100;
  char Buf[Buf_Length];                          // Buffer to store text from file
  unsigned uBuf;
  double lfBuf;
  Materials::Material Mat;

  // Now let's print the number of particles
  fprintf(File,   "       -- Particles --\n");
  fprintf(File,   "Number of particles:          %u\n\n",    Body_In.Get_Num_Particles());

  // We already have the Body's name, cuboid/boundary flags, and
  // dimensions (if cuboid). We can therefore skip over these lines.
  fgets(Buf, 99, File);                          // Skip 'name' line.
  fgets(Buf, 99, File);                          // Skip blank line
  fgets(Buf, 99, File);                          // Skip 'Is a cuboid' line

  if(Body_In.Get_Cuboid() == true) {           // if a cuboid, skip the 3 dimension lines.
    fgets(Buf, 99, File);
    fgets(Buf, 99, File);
    fgets(Buf, 99, File);
  } // if(Body_In.Get_Cuboid() == true) {

  fgets(Buf, 99, File);                          // Skip 'is a boundary' line
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

  // Set first time step flag to false
  Body_In.First_Time_Step == false;

  // Now read in number of particles
  fgets(Buf, 99, File);                          // Skip 'Particles' line.
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Num_Particles);

  // Now read in particles.
  for(unsigned i = 0; i < Num_Particles; i++)
    Load_Particle(Body_In[i], File, Body_In.Get_Cuboid());

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

  return 0;
} // int Data_Dump::Load_Data_From_File(unsigned & Num_Particles, Body Body_In) {



void Data_Dump::Load_Particle(Particle & P_In, FILE * File, const bool Is_Cuboid) {
  /* This function reads in the particle data for a specific particle. This
  data is transfered from the File to P_In. */

  char Buf[100];                                 // Buffer to store text from file (discarded)

  // First, read in the particle's ID and dimensions
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &P_In.ID);

  if(Is_Cuboid == true) {
    fread(Buf, 1, 30, File); fscanf(File, " %u %u %u\n", &P_In.i, &P_In.j, &P_In.k);
  }

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
  for(unsigned i = 0; i < P_In.Num_Neighbors; i++)
    fscanf(File, " %u", &P_In.Neighbor_IDs[i]);
} // void Data_Dump::Load_Particle(Particle & P_In, const FILE * File, const bool Is_Cuboid) {

#endif
