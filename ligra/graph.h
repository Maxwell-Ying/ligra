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

struct Empty_Mem : public Deletable {
  void del() {}
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

  vertex * getvertex(uintT j) {
    return &V[j];
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

  int get_version() {
    return version;
  }

  void update_m(void) {
    int ret = 0;
    for (auto v : V) {
      ret += v.outNeighbors.size();
    }
    m = ret;
  }

  void set_version(int vers) {
    version = vers;
  }

  uintE accessAllEdges() {
    uintE ret = 0;
    // int count = 0;
    // auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < n; i++) {
      for (auto j = 0; j < V[i].getInDegree(); j++) {
        ret ^= V[i].getInNeighbor(j);
      }
      // count += V[i].getInDegree();
      // if (!(i % 1000000) && i) {
        // auto finish = std::chrono::high_resolution_clock::now();
        // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";
        // cout << "count : " << count << endl;
      // }
    }
    // abort();
    return ret;
  }
};

template <>
struct graph<hybridVertex> {
  myVector<hybridVertex*> V;
  int version;
  long m;
  long n;
  bool transposed;
  uintE * flags;
  Deletable * D;

  graph(myVector<hybridVertex *> _V, long _n, long _m, Deletable* _D) : V(_V), n(_n), m(_m),
  D(_D), flags(NULL), transposed(0), version(0) {}

  void del() {
    if (flags != NULL) free(flags);
    D->del();
    free(D);
    if (V.size() > 0) {
      for (auto i = 0; i < V.size(); i++) {
        if (i) {
          delete V[i];
          V[i] = NULL;
        }
      }
    }
  }

  hybridVertex * getvertex() {
    return V[0];
  }

  hybridVertex * getvertex(uintT j) {
    if (j >= n) {
      cout << "wrong vertex index" << endl;
      abort();
    }
    return V[j];
  }

  int get_edge_number() {
    int ret = 0;
    for (auto i : V) {
      ret += i->getInDegree();
    }
    return ret;
  }
  
  int get_version() {
    return version;
  }

  void set_version(int vers) {
    version = vers;
  }
  uintE accessAllEdges() {
    uintE ret = 0;
    for (auto i: V) {
      for (auto j=0; j<i->getInDegree(); j++) {
        ret &= i->getInNeighbor(j);
      }
    }
    return ret;
  }
};
#endif
