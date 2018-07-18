#if !defined(VTK_FILE_SOURCE)
#define VTK_FILE_SOURCE

#include "VTK_File.h"
#include "Particle.h"

void VTK_File::Get_File_Name(string & Str) {
  char Buf[6];
  sprintf(Buf,"%05d",File_Number);
  File_Number++;

  Str += "_variables_";
  Str += Buf;
  Str += ".vtk";
} // void Get_File_Name(string & Str) {

void VTK_File::Add_Point_Data(FILE * File, char * Weight_Name, unsigned int Num_Particles, double * Data) {
  // Print header.
  fprintf(File,"SCALARS ");
  fprintf(File,Weight_Name);
  fprintf(File," float\n");
  fprintf(File,"LOOKUP_TABLE default\n");

  // Now print supplied data to file
  for(unsigned int i = 0; i < Num_Particles; i++) {
    fprintf(File,"\t %8.3f\n",Data[i]);
  } // for(unsigned int i = 0; i < Num_Particles; i++) {
}

void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles) {
  // Set up file
  string File_Name = "Test";
  Get_File_Name(File_Name);

  string File_Path = "./Position_Files/";
  File_Path += File_Name;

  FILE * File = fopen(File_Path.c_str(), "w");

  //////////////////////////////////////////////////////////////////////////////
  // Print file header
  fprintf(File,"%s\n","# vtk DataFile Version 3.0");
  fprintf(File,"%s\n","test_file");
  fprintf(File,"%s\n","ASCII");
  fprintf(File,"%s\n","DATASET POLYDATA");
  fprintf(File,"POINTS %i float\n",Num_Particles);

  //////////////////////////////////////////////////////////////////////////////
  // Cycle through particles, print spacial positions of each particle
  Vector x;
  for(unsigned int i = 0; i < Num_Particles; i++) {
    x = Particles[i].x;

    fprintf(File,"%8.3f \t %8.3f \t %8.3f\n",x[0], x[1], x[2]);
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  //////////////////////////////////////////////////////////////////////////////
  /* Find the components of S and E for each particle
  We Calculate each of these by first finding P and F for each particle. We then
  use these quantities to calculate Sigma and E. The components of each of these
  tensors are then stored in dynamic arrays. Once this is finished, we write the
  components to the output file.

  We choose to store the components in dynamic arrays before writing to the file
  to improve performnace. Each particle's P, F tensors are in distinct memory
  locations. If we were to collect each component at a time, we'd have to read
  in a new cache line for each component from each particle. In other words,
  we'd only get one double from each cache line. This is bad. To improve this,
  we only pull from the particles once. We pull in the entire tensor, then write
  its components to the dynamic arrays. We use all 6 of each tensor's components
  in each iteration. This means that we can make better usage of cache lines.

  Once the components have been written to the arrays, tensor components
  belonging to adjacent particles are next to each ther in memory. For example,
  in the P11 array, P11 component of the 100th particle is stored right next to
  the P11 component of the 101th particle. when we write to the file, we pull
  from the dynamic arrays. This allows us to use cache lines efficiently.

  This improves performance.
  */

  /* Create dynamic arrays for components of S, E (note, both are symmetric, so
  we only need to store 6 components) and J (det F)*/

  double * LamM = new double[Num_Particles];
  //double * LamH = new double[Num_Particles];
  //double * LamC = new double[Num_Particles];
  double * D = new double[Num_Particles];
  /*
  double * S11 = new double[Num_Particles];
  double * S22 = new double[Num_Particles];
  double * S33 = new double[Num_Particles];
  double * S21 = new double[Num_Particles];
  double * S31 = new double[Num_Particles];
  double * S32 = new double[Num_Particles];
  */
  //double * E11 = new double[Num_Particles];
  //double * E22 = new double[Num_Particles];
  //double * E33 = new double[Num_Particles];
  //double * E21 = new double[Num_Particles];
  //double * E31 = new double[Num_Particles];
  //double * E32 = new double[Num_Particles];

  //double * J = new double[Num_Particles];

  Tensor F, P, S, E;
  Tensor I{1,0,0,
           0,1,0,
           0,0,1};

  for(unsigned int i = 0; i < Num_Particles; i++) {
    LamM[i] = Particles[i].Stretch_M;
    //LamH[i] = Particles[i].Stretch_H;
    //LamC[i] = Particles[i].Stretch_Critical;
    D[i] = Particles[i].D;

    // Get F, P from current particle
    //F = Particles[i].F;
    //P = Particles[i].P;

    // Use F to calculate determinant (J)
    //J[i] = Determinant(F);

    // Calculate S from P.
    //S = P*(F^(T))/J[i];

    // Get components of S
    /*
    S11[i] = S[3*0 + 0];
    S22[i] = S[3*1 + 1];
    S33[i] = S[3*2 + 2];
    S21[i] = S[3*1 + 0];
    S31[i] = S[3*2 + 0];
    S32[i] = S[3*2 + 1];
    */

    // Now calculate E = (1/2)(C-I)
    //E = (1./2.)*((F^T)*F - I);

    // Now get components of E
    //E11[i] = E[3*0 + 0];
    //E22[i] = E[3*1 + 1];
    //E33[i] = E[3*2 + 2];
    //E21[i] = E[3*1 + 0];
    //E31[i] = E[3*2 + 0];
    //E32[i] = E[3*2 + 1];
  } // for(unsinged int i = 0; i < Num_Particles; i++) {

  // Now print these values to the file.
  fprintf(File,"POINT_DATA %i\n", Num_Particles);
  char Weight_Name[5];

  /* Damage paramaters */

  std::strcpy(Weight_Name, "LamM");
  Add_Point_Data(File, Weight_Name, Num_Particles, LamM);

  //std::strcpy(Weight_Name, "LamH");
  //Add_Point_Data(File, Weight_Name, Num_Particles, LamH);

  //std::strcpy(Weight_Name, "LamC");
  //Add_Point_Data(File, Weight_Name, Num_Particles, LamC);

  std::strcpy(Weight_Name, "D");
  Add_Point_Data(File, Weight_Name, Num_Particles, D);
  /* Components of S */
  /*
  std::strcpy(Weight_Name, "S11");
  Add_Point_Data(File, Weight_Name, Num_Particles, S11);

  std::strcpy(Weight_Name, "S22");
  Add_Point_Data(File, Weight_Name, Num_Particles, S22);

  std::strcpy(Weight_Name, "S33");
  Add_Point_Data(File, Weight_Name, Num_Particles, S33);

  std::strcpy(Weight_Name, "S21");
  Add_Point_Data(File, Weight_Name, Num_Particles, S21);

  std::strcpy(Weight_Name, "S31");
  Add_Point_Data(File, Weight_Name, Num_Particles, S31);

  std::strcpy(Weight_Name, "S32");
  Add_Point_Data(File, Weight_Name, Num_Particles, S32);
  */

  /* Components of E*/
  //std::strcpy(Weight_Name, "E11");
  //Add_Point_Data(File, Weight_Name, Num_Particles, E11);

  //std::strcpy(Weight_Name, "E22");
  //Add_Point_Data(File, Weight_Name, Num_Particles, E22);
  /*
  std::strcpy(Weight_Name, "E33");
  Add_Point_Data(File, Weight_Name, Num_Particles, E33);

  std::strcpy(Weight_Name, "E21");
  Add_Point_Data(File, Weight_Name, Num_Particles, E21);

  std::strcpy(Weight_Name, "E31");
  Add_Point_Data(File, Weight_Name, Num_Particles, E31);

  std::strcpy(Weight_Name, "E32");
  Add_Point_Data(File, Weight_Name, Num_Particles, E32);
  */

  /* J */
  /*
  std::strcpy(Weight_Name, "J");
  Add_Point_Data(File, Weight_Name, Num_Particles, J);
  */

  // Deallocate dynamic arrays
  delete [] LamM;
  //delete [] LamH;
  //delete [] LamC;
  delete [] D;
  /*
  delete [] S11;
  delete [] S22;
  delete [] S33;
  delete [] S21;
  delete [] S31;
  delete [] S32;
  */
  //delete [] E11;
  //delete [] E22;
  //delete [] E33;
  //delete [] E21;
  //delete [] E31;
  //delete [] E32;

  //delete [] J;

  // Free the file
  fclose(File);
} // void Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles) {

#endif
