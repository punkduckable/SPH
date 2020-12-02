#include "FEB_File.h"
#include "Vector/Vector.h"
#include "Errors.h"

void IO::Read_FEB_File(const std::string & File_Name, Vector ** X_Ptr, unsigned & Num_Nodes) {
  /* This function, as the name suggests, is used to read in an FEB file. It
  reads in the position of the nodes in that body and then stores them in the
  X_Ptr pointer. */

  // First, we need to get the path to the Febio file. Now we can open the file
  std::string File_Path = "./IO/Files/";
  File_Path += File_Name;                        // Add file name
  File_Path += ".feb";                           // add FEBio file extension

  // Now open up the file.
  FILE * File = fopen(File_Path.c_str(), "r");

  // Check that file could be opened
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Can't Open File Exception: Thrown by Read_FEB_File\n"
            "For some reason, ./IO/Files/%s.feb can not be opened",
            File_Name.c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {

  // Now search through FEBio file until we come to 'Nodes' section
  unsigned i;
  bool Nodes_Found = false;
  long Start_Of_Nodes;
  const unsigned Buf_Length = 100;
  char Buf[Buf_Length+1];

  // Assign null character to final element of Buf.
  Buf[Buf_Length] = '\0';

  while(fgets(Buf, Buf_Length, File) != nullptr && Nodes_Found == false) {
    for(i = 0; i < Buf_Length - 5; i++) {
      if(Buf[i] ==   '<' &&
         Buf[i+1] == 'N' &&
         Buf[i+2] == 'o' &&
         Buf[i+3] == 'd' &&
         Buf[i+4] == 'e' &&
         Buf[i+5] == 's') {

         Nodes_Found = true;
         break;
      } // if (Buf = 'Nodes')

      else if(Buf[i] == '>')
        break;
    } // for(i = 0; i < Buf_Length - 5; i++) {

    if(Nodes_Found == true) { break; }
  } // while(fgets(Buf, Buf_Length, File) != nullptr) {

  // Record current file position.
  Start_Of_Nodes = ftell(File);

  long Start_Of_Current_Line, End_Of_Current_Line;
  bool End_Of_Nodes = false;
  Num_Nodes = 0;

  // Now cycle through lines of code, reading in number of IDs
  while(End_Of_Nodes == false) {
    Start_Of_Current_Line = ftell(File);

    // Read in a line and check if we're at the end of the nodes section
    if(fgets(Buf, Buf_Length, File) == nullptr) {
      throw Bad_Read("Bad Read Exception: Thrown by IO::Read_FEB_File\n"
                     "The program could not find the end of the FEBio file :(\n");
    } // if(fgets(Buf, Buf_Length, File) == nullptr) {

    End_Of_Current_Line = ftell(File);

    for(i = 0; i < Buf_Length-7; i++) {
      if(Buf[i]   == '<' &&
         Buf[i+1] == '/' &&
         Buf[i+2] == 'N' &&
         Buf[i+3] == 'o' &&
         Buf[i+4] == 'd' &&
         Buf[i+5] == 'e' &&
         Buf[i+6] == 's' &&
         Buf[i+7] == '>') {

        End_Of_Nodes = true;
        break;
      } // if( Buf = <\Nodes>)

      if(Buf[i] == '\n') { break; }
    } // for(i = 0; i < Buf_Length) {

    // Exit for loop if we're done.
    if(End_Of_Nodes == true) { break; }

    // Read in Node ID.
    fseek(File, Start_Of_Current_Line, SEEK_SET);     // Move file pointer to start of line
    fscanf(File, " <node id=\"%u\">", &Num_Nodes);    // Read in ID
    fseek(File, End_Of_Current_Line, SEEK_SET);       // Move file pointer back to end of line
  } // while(End_Of_Nodes = false) {



  //////////////////////////////////////////////////////////////////////////////
  // Read in node positions.

  // Now that we know number of nodes, go back to start of nodes
  fseek(File, Start_Of_Nodes, SEEK_SET);

  // Allocate an array to hold the reference positons of the nodes
  unsigned ID_Buf;
  *X_Ptr = new Vector[Num_Nodes];

  for(i = 0; i < Num_Nodes; i++) {
    fscanf(File, " <node id=\"%u\">%le, %le, %le</node>\n", &ID_Buf, &(*X_Ptr)[i](0), &(*X_Ptr)[i](1), &(*X_Ptr)[i](2));
  } // for(i = 0; i < Num_Nodes; i++) {

  printf("Read in %u particles \n",ID_Buf);
  fclose(File);
} // void IO::Read_FEB_File(const std::string & File_Name, Vector ** X_Ptr, unsigned & Num_Nodes) {
