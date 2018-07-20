#if !defined(FEB_SOURCE)
#define FEB_SOURCE

int FEB_File::Read_FEB_File(const string & File_Name) {
  // First, we need to get the path to the Febio file. Now we can open the file
  string File_Path = "../Files/";
  File_Path += File_Name;

  // Now open up the file.
  FILE * File = fopen(File_Path.c_str(), "r");

  // Check that file could be opened
  if(File == NULL) {
    printf("Requested FEBio file does not exist.\n");
    return 1;
  } // if(File == NULL) {

  // Now search through FEBio file until we come to 'Nodes' section
  unsigned int i;
  bool Nodes_Found = false;
  long Start_Of_Nodes;
  const unsigned int Buf_Length = 100;
  char Buf[Buf_Length];

  while(fgets(Buf, Buf_Length, File) != NULL && Nodes_Found == false) {
    for(i = 0; i < Buf_Length - 5; i++) {
      if(Buf[i] == '<') {
        if(Buf[i+1] == 'N' &&
           Buf[i+2] == 'o' &&
           Buf[i+3] == 'd' &&
           Buf[i+4] == 'e' &&
           Buf[i+5] == 's') {

           Nodes_Found = true;
           break;
         } // if (Buf = 'Nodes')
      } // if(Buf[i] == '<') {

      else if(Buf[i] == '>')
        break;
    } // for(i = 0; i < Buf_Length - 5; i++) {
  } // while(fgets(Buf, Buf_Length, File) != NULL) {

  // Record current file position.
  Start_Of_Nodes = ftell(File);

  long Start_Of_Current_Line;
  bool End_Of_Nodes = false;
  unsigned int Num_Nodes;

  // Now cycle through lines of code, reading in number of IDs
  while(End_Of_Nodes == false) {
    Start_Of_Current_Line = ftell(File);

    // First check if we're at the end of the nodes section
    if(fgets(Buf, Buf_Length, File) == NULL) {
      printf("Couldn't find end of FEBio file\n");
      return 1;
    }

    for(i = 0; i < Buf_Length; i++) {
      if(Buf[i] == '<') {
        if(Buf[i+1] == '\\' &&
           Buf[i+2] == 'N') {

            End_Of_Nodes = true;
            break;
         } // if( Buf = <\Nodes>)
      } // if(Buf[i] == '<') {
    } // for(i = 0; i < Buf_Length) {

    // Exit for loop if we're done.
    if(End_Of_Nodes == true)
      break;

    // Read in Node ID.
    fseek(File, 0, Start_Of_Current_Line);
    fscanf(File, " <node id=\"%u\">", &Num_Nodes);
  } // while(End_Of_Nodes = false) {

  printf("Number of nodes = %u\n",Num_Nodes);

  // Now that we know number of nodes, go back to start of nodes
  fseek(File, 0, Start_Of_Nodes);

  return 0;
} // int FEB_File::Read_FEB_File(const string & File_Name) {

#endif
