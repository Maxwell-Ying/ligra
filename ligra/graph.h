#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>
#include "vertex.h"
#include "compressedVertex.h"
#include "parallel.h"
#include "myVector.h"
#include "get_mem.h"

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
private:
  int max_version;
  myVector<vertex> V;
  int version;//version id
public:
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

  graph(myVector<vertex>& _V, long _n, long _m, Deletable* _D) : V(_V), n(_n), m(_m),
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

  inline vertex * getvertex(uintT j) {
    if (j >= n) {
      resize(j+1);
    }
    return &V[j];
  }

  void resize(uintT new_size) {
    cout << "vertex expand from " << n << " to " << new_size << endl;
    V.resize(new_size);
    for (auto i=n; i<V.size(); i++) {
        V[i].prepare();
      }
    n = V.size();
  }

  int get_edge_number() {
    return m;
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

  int get_max_version() {
    return max_version;
  }

  void update_m(void) {
    int ret = 0;
    for (auto v : V) {
      ret += v.getInDegree();
    }
    m = ret;
  }

  void add_m(int a) {
    if (m+a < 0) {
      update_m();
    }
    m += a;
  }

  void set_version(int vers) {
    version = vers;
    if (version > max_version) {
      max_version = version;
    }
  }

  uintE accessAllEdges() {
    uintE ret = 0;
    int count = 0;
    int length = 0;
    for (auto i = 0; i < n; i++) {
      length = V[i].getInDegree();
      vertex & tmp = V[i];
      for (auto j = 0; j < length; j++) {
        ret ^= tmp.getInNeighbor(j);
        // auto k = tmp.getInNeighbor(j);
        // tmp.empty_step();
      }
      count += length;
    }
    cout << "count " << count << " edges " << endl;
    return ret;
  }

  uintE access_vertex(uintE vtx) {
    uintT ret = 0;
    for (auto i=0; i<V[vtx].getInDegree(); i++) {
      ret ^= V[vtx].getInNeighbor(i);
    }
    return ret;
  }

  myVector<uintE> get_ldvertex(uintT ver) {
    myVector<uintE> ret;
    for (auto i=0; i<n; i++) {
      if (V[i].is_ld_to(ver)) {
        ret.push_back(i);
      }
    }
    return ret;
  }

  bool is_hybrid() {
    return sizeof(vertex) == sizeof(HVertex);
  }

  bool is_high_degree_vertex(uintT index) {
    return V[index].is_high_degree_vertex();
  }
};

#endif
