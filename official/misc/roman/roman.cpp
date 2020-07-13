// Roman Numerals
//
// Author: Daniel Anderson
// Date: 18-01-2017
// Reliability: 5
// Tested on:
//   ECNA07_A, SPOJ_ROMANN, SPOJ_ROMAN008, UVA_759, UVA_11616,
//   UVA_12397
//

//listings:roman
// Arabic / Roman numeral conversion for 0 < x < 4000. Just be greedy from high to low.
const string R[13] = {"M","CM","D","CD","C","XC","L","XL","X","IX","V","IV","I"};
const int A[13] = {1000,900,500,400,100,90,50,40,10,9,5,4,1};

string to_roman(int x) {  // Convert x to Roman numerals (0 < x < 4000)
  string s;               // For x >= 4000, additional 'M's are appended
  for (int i=0; i<13; i++)
    while (x >= A[i])
      x -= A[i], s += R[i];
  return s;
}

int to_arabic(string s) {  // Convert the Roman numeral s into Arabic
  int x = 0;               // Additional leading 'M's will be treated as thousands.
  for (int i=0; i<13; i++)
    while (s.find(R[i])==0)
      x += A[i], s.erase(0,R[i].size());
  return x;
}
//listings:/roman
