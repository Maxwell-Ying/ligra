#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "vertex.h"
#include "compressedVertex.h"
#include "parallel.h"
#include "myVector.h"

using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

// Class that handles implementation specific freeing of memory
// owned by the graph
struct Deletable {
public:
  virtual void del() = 0;
};

template <class vertex>
struct Uncompressed_Mem : public Deletable {
public:
  myVector<vertex> V;
  long n;
  long m;
  void* allocatedInplace, * inEdges;

  Uncompressed_Mem(myVector<vertex> VV, long nn, long mm, void* ai, void* _inEdges = NULL)
  : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges) { }
  Uncompressed_Mem(vertex* VV, long nn, long mm, void* ai, void* _inEdges = NULL)
  :  n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges) {
    for(int i=0; i<nn; i++) {
      V.push_back(*(VV+i));
    }
  }
  void del() {
    if (allocatedInplace == NULL)
      for (long i=0; i < n; i++) V[i].del();
    // else free(allocatedInplace);
    V.clear();
    if(inEdges != NULL) free(inEdges);
  }
};

template <class vertex>
struct Compressed_Mem : public Deletable {
public:
  vertex* V;
  char* s;

  Compressed_Mem(vertex* _V, char* _s) :
                 V(_V), s(_s) { }

  void del() {
    free(V);
    free(s);
  }
};

template <class vertex>
struct graph {
  myVector<vertex> V;
  int version;//version id
  long n;
  long m;
  bool transposed;
  uintE* flags;
  Deletable *D;

graph(vertex* _V, long _n, long _m, Deletable* _D) : n(_n), m(_m),
  D(_D), flags(NULL), transposed(0), version(0) {
    for(int i=0; i<_n; i++) {
      V.push_back(*(_V+i));
    }
  }
graph(myVector<vertex> _V, long _n, long _m, Deletable* _D) : V(_V), n(_n), m(_m),
  D(_D), flags(NULL), transposed(0), version(0) {}

graph(vertex* _V, long _n, long _m, Deletable* _D, uintE* _flags, int _version) : V(_V),
  n(_n), m(_m), D(_D), flags(_flags), transposed(0), version(_version) {}

  void del() {
    if (flags != NULL) free(flags);
    D->del();
    free(D);
  }

  void transpose() {
    if ((sizeof(vertex) == sizeof(asymmetricVertex)) ||
        (sizeof(vertex) == sizeof(compressedAsymmetricVertex))) {
      parallel_for(long i=0;i<n;i++) {
        V[i].flipEdges();
      }
      transposed = !transposed;
    }
  }
  vertex* getvertex() {
	  return V.data();
  }

  int get_edge_number() {
    int count = 0;
    for (auto v : V) {
      count += v.outNeighbors.size();
    }
    return count;
  }

  int get_edge_capicity() {
    int count = 0;
    for (auto v : V) {
      count += v.outNeighbors.get_cap();
    }
    return count;
  }

  int getversion() {
    return version;
  }

  void setversion(int vers) {
    version = vers;
  }
};
#endif
