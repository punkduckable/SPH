#if !defined(VTK_FILE_SOURCE)
#define VTK_FILE_SOURCE

#include "VTK_File.h"
#include "Particle.h"

bool VTK_File::Append_Digits(unsigned int N, string & Str) {
  /* To convert N into a string, we need to get the digits of N.

  To do this, we first need to figure out the number of digits in N. Once we
  know that, we can dynamically allocate an array to store the digits of N.
  To get the actaul digits, we can get the digits of N using integer division!

  To get the 1's place of an int, we can use the following: N - 10*(N/10)
  (or just N % 10). Once we have that, we can divide N by 10 and repeat. This
  process continues until the remaining number is equal to zero. */

  // Get number of digits in N
  unsigned int Num_Digits = 0;
  unsigned int N_temp = N;

  do {
    Num_Digits++;
    N_temp /= 10;
  } while(N_temp != 0);

  // Check that we are less than the max number of digits.
  if(Num_Digits > File_Number_Max_Digits) {
    printf("File Counter is too high! \n");
    return true;
  }

  /* Now dynamically allocate a char array to store the digits of N. Note that
  we allocate an extra element in this string for the 'end-of-string'
  character */
  char * Digits = new char[File_Number_Max_Digits+1];
  Digits[File_Number_Max_Digits] = '\0';

  // Populate the elements of the Digits array
  for(int i = File_Number_Max_Digits-1; i >= File_Number_Max_Digits - Num_Digits; i--) {
    Digits[i] = N%10+48;
    N /= 10;
  } // for(int i = Num_Digits-1; i >= 0; i--) {

  for(int i = 0; i < File_Number_Max_Digits - Num_Digits; i++) {
    Digits[i] = '0';
  } // for(int i = 0; i < File_Number_Max_Digits - Num_Digits; i++) {

  // Now return the string version of N
  Str += Digits;
  return false;
} // bool Append_Digits(unsigned int N, string & Str) {

void VTK_File::Get_File_Name(string & Str) {
  File_Number++;
  Str += "_variables_";
  Append_Digits(File_Number, Str);
  Str += ".vtk";
} // void Get_File_Name(string & Str) {

void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles) {
  string File_Name = "Test";
  Get_File_Name(File_Name);

  string File_Path = "./Position_Files/";
  File_Path += File_Name;

  FILE * File = fopen(File_Path.c_str(), "w");

  // Print file header
  fprintf(File,"%s\n","# vtk DataFile Version 3.0");
  fprintf(File,"%s\n","test_file");
  fprintf(File,"%s\n","ASCII");
  fprintf(File,"%s\n","DATASET POLYDATA");
  fprintf(File,"POINTS %i float\n",Num_Particles);

  // Cycle through particles, print spacial positions of each particle
  Vector x;
  for(int i = 0; i < Num_Particles; i++) {
    x = Particles[i].Get_x();

    fprintf(File,"%10.5f, %10.5f, %10.5f;\n",x(0), x(1), x(2));
  }

  // Now print weighting information
  fprintf(File,"POINT_DATA %i\n", Num_Particles);
  fprintf(File,"SCALARS C float\n");
  fprintf(File,"LOOKUP_TABLE default\n",Num_Particles);

  double P_1_1;
  for(int i = 0; i < Num_Particles; i++) {
    P_1_1 = Particles[i].Get_P_1_1();
    fprintf(File,"\t %14.7f\n",P_1_1);
  } // for(int i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles) {

#endif
