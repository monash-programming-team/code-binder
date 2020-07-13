#include "test-includes.h"
#include "../official/misc/roman/roman.cpp"

TEST_CASE ("Brute force all Roman Numerals", "[roman]" ){
  for(int i=1;i<4000;i++){
    string roman = to_roman(i);
    REQUIRE( i == to_arabic(roman) );
  }
}
