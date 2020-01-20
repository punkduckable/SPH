#include "IO_OPS.h"
#include "Simulation/Simulation.h"
#include "Vector/Vector.h"
#include "Errors.h"

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
