#include "Array.h"
#include "List.h"

TEST_CASE("Array Tests","[Array]") {
  // First, make an array using the default constructor
  const unsigned Length = 7;
  Array<unsigned> A1{};

  REQUIRE( A1.Get_Length()     == 0);
  REQUIRE( A1.Get_Length_Set() == false);
  REQUIRE_THROWS( A1[1] );

  // Now set the length, test that it set up properly.
  A1.Set_Length(Length);
  REQUIRE( A1.Get_Length()     == Length);
  REQUIRE( A1.Get_Length_Set() == true);

  // Try setting the length again (this should throw an exception)
  REQUIRE_THROWS( A1.Set_Length(Length) );

  // Try setting some elements and accessing their values
  for(unsigned i = 0; i < Length; i++) { A1[i] = i; }
  for(unsigned i = 0; i < Length; i++) { REQUIRE( A1[i] == i); }



  // Make another list using the defualt constructor but is set up using Setup_From_List
  Array<unsigned> A2{};
  List<unsigned> L1;
  for(unsigned i = 0; i < Length; i++) { L1.Push_Back(i); }

  // Set up A2 from list
  A2.Setup_From_List(L1);
  REQUIRE( A2.Get_Length()     == Length);
  REQUIRE( A2.Get_Length_Set() == true);

  // Make sure that A2 was set correctly.
  for(unsigned i = 0;  i < Length; i++) { REQUIRE( A2[i] == i); }

  // Try setting up the array again (this should throw an exception)
  REQUIRE_THROWS( A2.Set_Length(Length) );
  REQUIRE_THROWS( A2.Setup_From_List(L1) );



  // Now, make an array with a preset number of elements.
  Array<unsigned> A3{Length};

  REQUIRE( A3.Get_Length()     == Length);
  REQUIRE( A3.Get_Length_Set() == true);

  // Try setting the array length (this should throw an error).
  REQUIRE_THROWS( A2.Set_Length(Length) );



  // Now make an Array from a list.
  List<unsigned> L2;
  for(unsigned i = 0; i < Length; i++) { L2.Push_Back(i); }
  Array<unsigned> A4{L2};

  // Make sure that A4 was set correctly.
  REQUIRE( A4.Get_Length()     == Length);
  REQUIRE( A4.Get_Length_Set() == true);
  for(unsigned i = 0;  i < Length; i++) { REQUIRE( A4[i] == i); }

  // Try setting up the array again (this should throw an exception)
  REQUIRE_THROWS( A2.Set_Length(Length) );
  REQUIRE_THROWS( A2.Setup_From_List(L2) );
} // TEST_CASE("Array Tests","[Array]") {
