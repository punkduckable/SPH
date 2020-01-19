#include "IO_OPS.h"
#include "Simulation.h"

std::string read_line_after_char(FILE * File, const char delim) {
  /* This function scans through the File until the character 'delim' is found
  (or the end of file is reached). Once the char has been found, this function
  reads the contents of the File until the end of the current line. */

  // First, read until we find the delim character (or reach the end of file)
  char c;
  do {
    c = fgetc(File);
  } while(c != delim && c != EOF);

  // Now read the rest until we encounter a newline character (or eof)
  std::string s;
  c = fgetc(File);
  while(c != EOF) {
    s += c;
    if(c == '\n') { break; }
    else { c = fgetc(File); }
  } // while(c != EOF) {

  return s;
} // std::string read_line_after_char(FILE * File, const char c) {



bool String_Ops::Contains(const char* Buffer, const char* Word, unsigned Start_At) {
  /* Function description:
  This function determines if Word is contained in Buffer. The optional Start_At
  argument can be used to only search through part of the string. If, for
  example, Start_At = 5, then this function will only search for matches that
  begin at index 5 (or later) of Buffer.
  My inp reader frequently checks if a particular word is in a string. Thus, I
  wrote this function to automate that process. */

  /* Assumptions:
  This function assumes that both Buffer and Word are NULL TERMINATED strings.
  That is, I assume that both end with the \0 character. */

  /* First, check if Buffer has fewer than Start_At characters (which happens
  if there is a \0 in a index whose value is less than Start_At). If not, return
  false */
  for(unsigned i = 0; i < Start_At; i++) {
    if(Buffer[i] == '\0') { return false; }
  } // for(unsigned i = 0; i < Start_At; i++) {

  // Loop through the characters of Buffer.
  unsigned i = Start_At;
  while(Buffer[i] != '\0') {
    // At each one, see if Word starts at that character.
    unsigned j = 0;
    while(Buffer[i+j] == Word[j]) {
      j++;

      /* If we're still in here and we've reached the end of "Word" then
      we've found a match! */
      if(Word[j] == '\0') { return true; }

      /* If we haven't reached the end of Word but we have reached the end of
      Buffer then Buffer does not contain Word. */
      if(Buffer[i+j] == '\0') { return false; }
    } // while(Buffer[i+j] == Word[j]) {

    i++;
  } // while(Buffer[i] != '\0') {

  /* If we get here then we cycled through Buffer without finding a match.
  Thus, buffer does not contain Word. */
  return false;
} // bool String_Ops::Contains(const char* Buffer, const char* Word, unsigned Start_At) {



std::vector<std::string> String_Ops::Split(std::string & S, const char delim) {
  /* Function description:
  This function is designed to read in a string and split it using a delimeter.
  This works using the char array version of Split (see below).
  Delim is a defaulted argument. By default, delim = ','. */

  return String_Ops::Split(S.c_str(), delim);
} // std::vector<std::string> String_Ops:Split(std::string & S, const char delim) {



std::vector<std::string> String_Ops::Split(const char* S, const char delim) {
  /* Function description:
  This function is designed to read in a nul terminated char array and split it
  into a set of substrings. The splits occur whenever a delim character is found
  The characters between any two instances of delim (as well as any
  characters before the first instance of delim and the characters after the
  last instance of delim) are packaged together as a substring and added to the
  substring vector (which is what get's returned).
  Delim is a defaulted argument. By default, delim = ','. */

  std::vector<std::string> Sub_Strings;

  unsigned Index_Start = 0;                      // Index of the start of the current substring.
  unsigned N_Chars_Since_Delim = 0;              // Number of characters after the start of the substring that are not Delim

  unsigned i = 0;
  while(S[i] != '\0' && S[i] != '\n' && S[i] != '\r') {
    if(S[i] == delim) {
      // Make the new substring.
      std::string Sub_Str;
      for(unsigned j = 0; j < N_Chars_Since_Delim; j++) { Sub_Str += S[Index_Start + j]; }
      Sub_Strings.push_back(Sub_Str);

      // Update Index_Start (to just 1 character after the delim index)
      Index_Start = i+1;

      // We're starting a new substring, so N_Chars_Since_Delim = 0;
      N_Chars_Since_Delim = 0;
    } // if(S[i] == delim) {
    else { N_Chars_Since_Delim++; }

    i++;
  } // while(S[i] != '\0' && S[i] != '\n' && S[i] != '\r') {

  /* Finally, make a substring from the charcters after the last instance of
  delim. Note that this only happens if Index_Start < i. (since, at this point,
  i is the index of the \0 character for the string) */
  if(Index_Start < i) {
    std::string Sub_Str;
    for(unsigned j = 0; j < N_Chars_Since_Delim; j++) { Sub_Str += S[Index_Start + j]; }
    Sub_Strings.push_back(Sub_Str);
  } // if(Index_Start < i) {

  return Sub_Strings;
} // std::vector<std::string> String_Ops:Split(const char* S, const char delim) {



Vector read_box_BC(std::string BC_Str) {
  /* The BC string should be of the form "<BC_x, BC_y, BC_z>" where each of
  C_x/y/z are either doubles or the keyword "Free". This function reads in
  those values and parses them as a boundary condition, which is then returned
  as a vector. */

  /* First, let's trim the string. By trim, I mean remove everything except for
  what's between < and > */

  unsigned num_chars = BC_Str.length();
  std::string trimmed_BC_Str;
  bool start = false;
  for(unsigned i = 0; i < num_chars; i++) {
    if(BC_Str[i] == '<') { start = true; }
    if(BC_Str[i] == '>') { break; }

    if(start == true) { trimed_BC_Str += BC_Str[i]; }
  } // for(unsigned i = 0; i < num_chars; i++) {

  /* Now, split the trimmed_Bc_Str at commas */
  std::vector<std::string> Split_BC_Str = String_Ops::split(trimed_BC_Str, ',');

  /* If there are not three strings in the Split_BC_Str then complain */
  if(Split_BC_Str.length() != 3) {
    printf("Bad BC Read. Thrown by read_box_BC\n");
    abort();
  } // if(Split_BC_Str.length() != 3) {

  Vector BC;
  for(unsigned i = 0; i < 3; i++) {
    if(Str_Ops::Contains(Split_BC_Str[i], "Free")) { BC[i] = Simulation::Free_BC_Box; }
    else { sscanf(Split_Bc_Str[i], " %lf ", &BC[i]); }
  } // for(unsigned i = 0; i < 3; i++) {

  return BC;
} // Vector read_box_BC(std::string BC_Str) {
