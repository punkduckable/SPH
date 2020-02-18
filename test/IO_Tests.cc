#include "IO/IO_Ops.h"
#include <fstream>
#include <string>
#include <stdio.h>

TEST_CASE("read_line_after","[IO]") {
  /* First, open the test file. */
  std::ifstream File;
  File.open("./IO/Test/test_file.txt");

  REQUIRE( File.is_open() == true );

  /* Now read in stuff from the test file.
  We should have
      Value 1: 1
      Value 2: 2
      Value 3: 3.0
      Value 4: 4321 */
  unsigned uBuf;
  float fBuf;
  std::string strBuf;



  //////////////////////////////////////////////////////////////////////////////
  // Try to recover values 1-4
  strBuf = IO::read_line_after(File, "Value 1:");
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  REQUIRE( uBuf == 1 );


  strBuf = IO::read_line_after(File, "Value 2:");
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  REQUIRE( uBuf == 2 );


  strBuf = IO::read_line_after(File, "Value 3:");
  sscanf(strBuf.c_str(), " %f \n", &fBuf);
  REQUIRE( fBuf == 3.0 );


  strBuf = IO::read_line_after(File, "Value 4:");
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  REQUIRE( uBuf == 4321 );



  //////////////////////////////////////////////////////////////////////////////
  // Try to recover value 5 (which does not exist)
  REQUIRE_THROWS(IO::read_line_after(File, "Value 5:") );



  //////////////////////////////////////////////////////////////////////////////
  // Try to recover values out of order

  File.close();
  File.open("./IO/Test/test_file.txt");

  strBuf = IO::read_line_after(File, "Value 2:");
  REQUIRE_THROWS( IO::read_line_after(File, "Value 1:") );



  //////////////////////////////////////////////////////////////////////////////
  // Try case-insensitive reading.

  File.close();
  File.open("./IO/Test/test_file.txt");

  strBuf = IO::read_line_after(File, "vAlUe 1:", false);
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  REQUIRE( uBuf == 1);

  strBuf = IO::read_line_after(File, "VALUE 2:", false);
  sscanf(strBuf.c_str(), " %u \n", &uBuf);
  REQUIRE( uBuf == 2 );

  REQUIRE_THROWS( IO::read_line_after(File, "vAlue 3:", true) );
} // TEST_CASE("read_line_after","[IO]") {
