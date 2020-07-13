#ifndef DEBUG_H
#define DEBUG_H

#pragma GCC diagnostic ignored "-Wunused-function"

// ---------------------------------------------------
// DEBUGGING output MACROS
// --------------------------------------------------
//Output a comma separated list of variable names and their contents
#define DEBUG(args...) cerr << "[Line " << __LINE__ << "]: "; \
  { vector<string> _v = __split(#args, ','); __ferr(_v.begin(), args); }

//Output the contents of a collection, one value per line
#define DEBUG_2D(A) cerr << "[Line " << __LINE__ << "]: " << #A << " = \n"; \
  for (const auto& R : (A)) { cerr << '\t' << R << '\n'; }

vector<string> __split(const string& s, char c) {
  vector<string> v;
  stringstream ss(s);
  string x;
  while (getline(ss, x, c))
    v.emplace_back(x);
  return v;
}

void __err(vector<string>::iterator it) { cerr << endl; }
template<typename T, typename... Args>
void __err(vector<string>::iterator it, T a, Args... args) {
  cerr << ", " << it -> substr((*it)[0] == ' ', it -> length()) << " = " << a;
  __err(++it, args...);
}
template<typename T, typename... Args>
void __ferr(vector<string>::iterator it, T a, Args... args) {
  cerr << it -> substr((*it)[0] == ' ', it -> length()) << " = " << a;
  __err(++it, args...);
}
  
// ---------------------------------------------------
// Stream overloads
// --------------------------------------------------
//Pair
template<typename U, typename V>
 ostream& operator<<(ostream &s, const pair<U, V> &x) {
  s << "(" << x.first << ", " << x.second << ")";
  return s;
}

//Vector
template<typename U>
ostream& operator<<(ostream &s, const vector<U> &x) {
  s << "[";
  bool was = false;
  for (auto it : x) {
    if (was) s << ", ";
    was = true;
    s << it;
  }
  s << "]";
  return s;
}

// Deque
template<typename U>
ostream& operator<<(ostream &s, const deque<U> &x) {
  s << "[";
  bool was = false;
  for (auto it : x) {
    if (was) s << ", ";
    was = true;
    s << it;
  }
  s << "]";
  return s;
}

//Map
template<typename U, typename V>
ostream& operator<<(ostream &s, const map<U, V> &x) {
  s << "{";
  bool was = false;
  for (auto it : x) {
    if (was) s << ", ";
    was = true;
    s << it;
  }
  s << "}";
  return s;
}

//Set
template<typename U>
ostream& operator<<(ostream &s, const set<U> &x) {
  s << "{";
  bool was = false;
  for (auto it : x) {
    if (was) s << ", ";
    was = true;
    s << it;
  }
  s << "}";
  return s;
}

//MultiSet
template<typename U>
ostream& operator<<(ostream &s, const multiset<U> &x) {
  s << "{";
  bool was = false;
  for (auto it : x) {
    if (was) s << ", ";
    was = true;
    s << it;
  }
  s << "}";
  return s;
}

#endif  // DEBUG_H
