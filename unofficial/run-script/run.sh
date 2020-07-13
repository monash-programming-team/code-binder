// setup

nano/pico run t.c

#!/bin/bash		// might not need it (test)
clear ; clear
g++ $1 -g -std=c++0x -Wall -Wconversion -Wshadow -o sol || exit

for i in *.in; do
   echo --- $i
   ./sol < $i > o && (diff -y o ${i::-3}.[ao]?? > t || cat t) || cat o
done

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;
typedef vector<int> vi;

int main()
{
	return 0;
}

#!/bin/bash		// can enter into terminal instead
chmod +x run
for i in {A..Z}; do
	mkdir $i
	//# cp sample/$i*/* $i # customize for use with work station
	cp t.c $i/$i.cc
	ln -st $i ../run
done
ls *
A/run A/A*

setxkbmap us,us -variant ,colemak grp:rwin_toggle

// setup gedit
Edit/Preferences
View tab:
	Display Line Numbers
	Highlight current line
	Highlight matching brackets
Editor tab:
	Tab width = 4
	Enable automatic indentation
