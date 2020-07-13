# MONASH UNIVERSITY CODE BINDER #

This is a collection of algorithms and data structures mainly used for
competitive programming competitions. Any Monash student may submit
code to this repository.

## GENERAL WARNING ##

Use files at your own risk! We cannot guarantee anything about their
correctness (but we hope that any code here is correct). You are not
allowed to use code from this binder for your coursework.

## DIRECTORIES ##

### "official" (Official Code Binder) ###

* This is Monash's Official Code Binder.
* All files must meet the guidelines below.
* Only the coaches may put files in this folder.

### "unofficial" (Unofficial Code Binder) ###

* The author of files in this folder feel that this code is worthy of being in the code binder.
* If the code is high enough quality, the coaches will move it to the Official Code Binder (see below for more details).
* All files must meet the guidelines below.
* Anyone may put files in this folder.
* Please only edit your own files (or have permission from the owner).

### "dev" (Under Development) ###

* This is for work in progress (Reliability Factors of 3 or less)
* Anyone can put files in this folder.
* Please only edit your own files (or have permission from the owner).
* Files in this folder do *not* have to meet the guidelines below.

## GUIDELINES ##

The file must

* Be compilable as is (with all necessary includes/imports).
* Include a main() which shows a simple example of how the code works.
* Give sufficient comments at the top to:

1. Describe the problem the algorithm is attempting to solve.
2. Give the complexity of the algorithm.
3. Explain any place that the code may need to be altered. The changes should be small and easy.
4. Give a Reliability Factor.

## RELIABILITY FACTOR ##

Each piece of code needs a reliability factor:

* 0 - New code, barely tested. Use at your own risk.
* 1 or 2 - AC on 1 or 2 problem. Use at your own risk.
* 3 or 4 - AC on 3 or 4 problems. Still a little risky to use.
* 5 - AC on 5+ problems.
* 6 - Coach’s approval (Official Code Binder).

Note that if you generate lots (millions or billions) of random tests and have a way to verify your answer, this will be considered AC on an additional problem.

## CODING STYLE ##

* Maximum number of characters per line: 80
* Do not deal with dynamic memory or bit-level memory allocations. This means you should not use: (unless you actually need them) memcpy, malloc, calloc, new, delete, free, pointers, etc.
* The important thing is that it must be usable by a wide range of people, not just you!
* Variable names of parameters must be useful and easy to understand. Variable names that are only used inside the function do not need to be useful.
* Limit the use of global variables. (Using a class helps!)

## COMMENTS IN CODE ##

At the top of your code, please include the following header
(obviously replaced with the description of your algorithm!).


```
#!c++

// Dijkstra's Algorithm
//
// Author      : Darcy Best
// Date        : August 24, 2014
// Reliability : 2
// Tested On   : 10986
//
// Computes the shortest distance between a source vertex and all
// other vertices in a weighted graph (directed or undirected).
// 
// All edge weights must be non-negative.
//
// Complexity: O( E log E )
```

Note:

* For the complexity, explain what each variable is if it is not clear.
* For graphs, V is the number of nodes (vertices), E is the number of edges.

Above your function, include a description of the input and output
parameters. For example:


```
#!c++

// Return:
//  pair of vectors. The first is D. The second is P.
//  D   : D[v] is the minimum distance from src to v. If you cannot
//        reach v, then D[v] is undefined. (D is for “Distance”)
//  P   : P[v] is the vertex you visit just before v on a shortest path
//        from src to v. P[src] == src if there is no negative cycle
//        reachable from the source. If you cannot reach v, then
//        P[v] == -1. (P is for “Predecessor”). P[src]==INT_MAX if there
//        is a negative cycle reachable from the source.
//
// Parameters:
//  G   : The graph
//  src : The source vertex
//
// Note: D and P will be resized inside the function.
pair<vector<int>, vector<int>> bellmanford(const Graph& G,int src){
  ( ... Code ... )
}
```

## HOW TO SUBMIT YOUR CODE ##

If you want to submit your code, fork the repo and send pull request.

Note

* *We will not accept code that has a low Reliability Factor.*
* We will only have one version of each algorithm in the official code binder.
* Some algorithms have multiple implementations and we will allow at most one of each implementation in these cases (for example, network flow can be solved using Ford-Fulkerson, Edmonds-Karp, Dinic, Push-Relabel).