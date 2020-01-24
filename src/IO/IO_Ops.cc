#include "IO_OPS.h"
#include "Simulation/Simulation.h"
#include "Vector/Vector.h"
#include "assert.h"
#include "Errors.h"

std::string IO::read_line_after_char(FILE * File, const char delim) {
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



std::string IO::read_line_after(std::ifstream & File, const char* Word) {
  /* Read in the file line by line until one is found that contains the
  word (or the end of file is reached) */
  assert(File.is_open() == true);


  unsigned Buffer_Size = 265;
  char Line_Buffer[Buffer_Size];
  int index;
  while(true) {
    File.getline(Line_Buffer, Buffer_Size);
    index = IO::String_Ops::Index_After_Word(Line_Buffer, Word);

    /* If Word is in Line then index will be some positive number. Otherwise,
    index = -1. */
    if(index != -1) { break; }

    /* Now check if the end of file has been reached. If so then the requested
    word could not be found in the File and an exception must be thrown. */
    if(File.eof() == true) {
      char Buffer[500];
      sprintf(Buffer,
              "Bad Read Exception: Thrown by IO::read_line_after\n"
              "You tried to find the word \"%s\" in a file.\n"
              "However, this word could not be found (or was located in a part \n"
              "of the file that you already read)\n",
              Word);
      throw Bad_Read(Buffer);
    } // if(File.eof() == true) {
  } // while(true) {

  /* If we are here then Line_Buffer stores a line of the File that contains
  the Word. Moreover, index stores the index of the first character in
  Line_Buffer after the end of the word. We can now construct the string (to
  return) from this */
  std::string Trimmed_Line;
  for(unsigned k = index; Line_Buffer[k] != '\0'; k++) { Trimmed_Line += Line_Buffer[k]; }
  return Trimmed_Line;
} // std::string IO::read_line_after(std::ifstream & File, const char* Word) {



bool IO::String_Ops::Contains(const char* Buffer, const char* Word, unsigned Start_At) {
  /* Function description:
  This function determines if Word is contained in Buffer. The optional Start_At
  argument can be used to only search through part of the string. If, for
  example, Start_At = 5, then this function will only search for matches that
  begin at index 5 (or later) of Buffer. */

  if(IO::String_Ops::Index_After_Word(Buffer,Word,Start_At) == -1) { return false; }
  else { return true; }
} // bool String_Ops::Contains(const char* Buffer, const char* Word, unsigned Start_At) {



int IO::String_Ops::Index_After_Word(const char* Buffer, const char* Word, unsigned Start_At) {
  /* Function description:
  This function attempts to find Word in Buffer. If Word is in Buffer, then this
  function will return the index of the first character in Buffer past the
  end of the Word. If Word is not contained in Buffer, then this function returns
  -1.

  The optional Start_At argument tells the function to only search through
  the string after the Start_At'th character. If, for example, Start_At = 5
  then this function will only search for Word in Buffer after the 5th
  character of Buffer. */


  /* Assumptions:
  This function assumes that both Buffer and Word are NULL TERMINATED strings.
  That is, I assume that both end with the \0 character. */

  /* First, check if Buffer has fewer than Start_At characters (which happens
  if there is a \0 in a index whose value is less than Start_At). If not, return
  -1 (Word can't be in Buffer) */
  for(unsigned i = 0; i < Start_At; i++) {
    if(Buffer[i] == '\0') { return -1; }
  } // for(unsigned i = 0; i < Start_At; i++) {

  // Loop through the characters of Buffer.
  unsigned i = Start_At;
  unsigned j;
  while(Buffer[i] != '\0') {
    // At each one, see if Word starts at that character.
    j = 0;
    while(Buffer[i+j] == Word[j]) {
      j++;

      /* If we're still in here and we've reached the end of "Word" then
      we've found a match! i+j is the index of the first character beyond the
      end of Word in Buffer. */
      if(Word[j] == '\0') { return i+j; }

      /* If we haven't reached the end of Word but we have reached the end of
      Buffer then Buffer does not contain Word. Return -1 */
      if(Buffer[i+j] == '\0') { return -1; }
    } // while(Buffer[i+j] == Word[j]) {
    i++;
  } // while(Buffer[i] != '\0') {

  /* If we get here then we cycled through Buffer without finding a match.
  Thus, buffer does not contain Word. Return -1 */
  return -1;
} // int IO::String_Ops::Index_After_Word(const char* Buffer, const char* Word, unsigned Start_At) {
