#include "List.h"


TEST_CASE("List Tests","[List]") {
  // First, create a list and test that the two ends point to nothing
  List<unsigned> L;
  REQUIRE_THROWS( L.Front() );
  REQUIRE_THROWS( L.Back() );
  REQUIRE( L.Get_Length()  == 0 );

  const unsigned Value1 = 7;
  unsigned Value2;

  // Try Push/Pop_Back
  L.Push_Back(Value1);
  REQUIRE( L.Front()         == Value1 );
  REQUIRE( L.Back()          == Value1 );
  REQUIRE( L.Get_Length() == 1 );

  Value2 = L.Pop_Back();
  REQUIRE( Value2            == Value1 );
  REQUIRE_THROWS( L.Front() );
  REQUIRE_THROWS( L.Back() );
  REQUIRE( L.Get_Length()  == 0 );

  // Try poping an element from the back of an empty list
  REQUIRE_THROWS( L.Pop_Back() );



  // Try Push/Pop_Front
  L.Push_Front(Value1);
  REQUIRE( L.Front()         == Value1 );
  REQUIRE( L.Back()          == Value1 );
  REQUIRE( L.Get_Length() == 1 );

  Value2 = 0;
  Value2 = L.Pop_Front();
  REQUIRE( Value2             == Value1 );
  REQUIRE( L.Get_Length()  == 0 );
  REQUIRE_THROWS( L.Front() );
  REQUIRE_THROWS( L.Back() );

  // Try poping an element from the front of an empty list
  REQUIRE_THROWS( L.Pop_Front() );



  // Now try pushing multiple elements
  L.Push_Back(3);
  L.Push_Front(2);
  L.Push_Back(4);
  L.Push_Front(1);
  L.Push_Back(5);
  L.Push_Front(0);

  for(unsigned i = 0; i <= 5; i++) {
    Value2 = L.Pop_Front();
    REQUIRE( L.Get_Length()  == 5-i );
    REQUIRE( Value2 == i );
  } // for(unsigned i = 0; i <= 5; i++) {
} // TEST_CASE("List Tests","[List]") {
