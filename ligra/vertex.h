#ifndef VERTEX_H
#define VERTEX_H

#include <chrono>
#include "vertexSubset.h"
#include "myVector.h"
using namespace std;

namespace decode_uncompressed {

  // Used by edgeMapDense. Callers ensure cond(v_id). For each vertex, decode
  // its in-edges, and check to see whether this neighbor is in the current
  // frontier, calling update if it is. If processing the edges sequentially,
  // break once !cond(v_id).
  template <class vertex, class F, class G, class VS>
  inline void decodeInNghBreakEarly(vertex* v, long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    uintE d = v->getInDegree();
    if (!parallel || d < 1000) {
      for (size_t j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
          auto m = f.update(ngh, v_id);
#else
          auto m = f.update(ngh, v_id, v->getInWeight(j));
#endif
          g(v_id, m);
        }
        if(!f.cond(v_id)) break;
      }
    } else {
      parallel_for(size_t j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
          auto m = f.updateAtomic(ngh, v_id);
#else
          auto m = f.updateAtomic(ngh, v_id, v->getInWeight(j));
#endif
          g(v_id, m);
        }
      }
    }
  }

  // Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNgh(V* v, long i, F &f, G &g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
      auto m = f.updateAtomic(i,ngh);
#else
      auto m = f.updateAtomic(i,ngh,v->getOutWeight(j));
#endif
        g(ngh, m);
      }
    });
  }

  // Used by edgeMapSparse. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNghSparse(V* v, long i, uintT o, F &f, G &g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
        auto m = f.updateAtomic(i, ngh);
#else
        auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
        g(ngh, o+j, m);
      } else {
        g(ngh, o+j);
      }
    });
  }

  // Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
  // and compactly write all neighbors satisfying g().
  template <class V, class F, class G>
  inline size_t decodeOutNghSparseSeq(V* v, long i, uintT o, F &f, G &g) {
    uintE d = v->getOutDegree();
    size_t k = 0;
    for (size_t j=0; j<d; j++) {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
#ifndef WEIGHTED
        auto m = f.updateAtomic(i, ngh);
#else
        auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    }
    return k;
  }

  // Decode the out-neighbors of v, and return the number of neighbors
  // that satisfy f.
  template <class V, class F>
  inline size_t countOutNgh(V* v, long vtx_id, F& f) {
    uintE d = v->getOutDegree();
    if (d < 2000) {
      size_t ct = 0;
      for (size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
#ifndef WEIGHTED
        if (f(vtx_id, ngh))
#else
        if (f(vtx_id, ngh, v->getOutWeight(i)))
#endif
          ct++;
      }
      return ct;
    } else {
      size_t b_size = 2000;
      size_t blocks = 1 + ((d-1)/b_size);
      auto cts = array_imap<uintE>(blocks, [&] (size_t i) { return 0; });
      parallel_for_1(size_t i=0; i<blocks; i++) {
        size_t s = b_size*i;
        size_t e = std::min(s + b_size, (size_t)d);
        uintE ct = 0;
        for (size_t j = s; j < e; j++) {
          uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
          if (f(vtx_id, ngh))
#else
          if (f(vtx_id, ngh, v->getOutNeighbor(j)))
#endif
            ct++;
        }
        cts[i] = ct;
      }
      size_t count = 0;
      return pbbs::reduce_add(cts);
    }
  }

  // Decode the out-neighbors of v. Apply f(src, ngh) and store the result
  // using g.
  template <class V, class E, class F, class G>
  inline void copyOutNgh(V* v, long src, uintT o, F& f, G& g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
#ifdef WEIGHTED
      E val = f(src, ngh, v->getOutWeight(j));
#else
      E val = f(src, ngh);
#endif
      g(ngh, o+j, val);
    });
  }

  // TODO(laxmand): Add support for weighted graphs.
  template <class V, class Pred>
  inline size_t packOutNgh(V* v, long vtx_id, Pred& p, bool* bits, uintE* tmp) {
    uintE d = v->getOutDegree();
    if (d < 5000) {
      size_t k = 0;
      for (size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        if (p(vtx_id, ngh)) {
          v->setOutNeighbor(k, ngh);
          k++;
        }
      }
      v->setOutDegree(k);
      return k;
    } else {
      parallel_for(size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        tmp[i] = ngh;
        bits[i] = p(vtx_id, ngh);
      }
      size_t k = sequence::pack(tmp, v->getOutNeighbors(), bits, d);
      v->setOutDegree(k);
      return k;
    }
  }

}

struct HVertex
{
  myVector<uintE> outNeighbors;
  myVector<myVector<uintE>> edges;
  int pos;
  uintT degree;

  HVertex() {
    pos = -2;
    degree = 0;
  }

  HVertex(myVector<uintE> n, int _pos=-1) {
    if (_pos == -1) {
      outNeighbors = n;
      pos = -1;
      degree = n.size();
    } else {
      myVector<uintE> new_ver(n);
      outNeighbors.push_back(0);
      pos = 0;
      edges.push_back(new_ver);
      degree = n.size();
    }
  }

  HVertex(uintE *data, uintT length, int _pos = -1) {
    if (_pos = -1) {
      pos = -1;
      outNeighbors.push_back(data, data + length);
    } else {
      myVector<uintE> new_ver(data, data+length);
      outNeighbors.push_back(0);
      pos = 0;
      edges.push_back(new_ver);
    }
    degree = length;
  }

  void setInNeighbors(uintE* _i) { return ; }
  void setOutNeighbors(uintE* _i) { return ; }
  void setInDegree(uintT _d) {   }
  void setOutDegree(uintT _d) { setInDegree(_d);  }
  void flipEdges() {}

  uintE* getInNeighbors() {
    if (pos != -1) {
      return edges[pos].data();
    } else {
      return outNeighbors.data();
    }
  }
  uintE * getOutNeighbors() { return getInNeighbors(); }

  myVector<uintE> getLastVersion() {
    if (pos >= 0) {
      return edges[edges.size()-1];
    } else if (edges.empty()){
      auto tmp = myVector<uintE>();
      return tmp;
    }
  }

  void index_delete(uintT index) {
    if (pos == -1) {
      outNeighbors.index_delete(index);
      degree --;
    }
  }
  void index_addtion(uintE data, uintT index) {
    if (pos == -1) {
      outNeighbors.index_addtion(data, index);
      degree ++;
    }
  }
  void push_back(uintE data) {
    if (pos == -1 || pos == -2) {
      outNeighbors.push_back(data);
      degree ++;
      if (pos == -2) {
        pos = -1;
      }
    }
  }
  void push_back(uintE * data, uintT length) {
    if (pos == -1 || pos == -2) {
      outNeighbors.push_back(data, data+length);
      degree += length;
      if (pos == -2) {
        pos = -1;
      }
    }
  }
  void push_back(myVector<uintE>::iterator start, myVector<uintE>::iterator end) {
    if (pos == -1 || pos == -2) {
      outNeighbors.push_back(start, end);
      degree = outNeighbors.size();
      if (pos == -2) {
        pos = -1;
      }
    }
  }

  void pop_back() {
    if (pos == -1) {
      outNeighbors.pop_back();
      degree --;
    }
  }

  bool is_high_degree_vertex() {
    return pos == -1;
  }

  void build(myVector<uintE> data) {
    if (pos != -2) {
      cout << "build from exist version " << endl;
    }
    outNeighbors.push_back(0);
    myVector<uintE> new_ver(data);
    edges.push_back(new_ver);
    pos = 0;
    degree = data.size();
  }

  void build(uintE* data, uintT length) {
    if (pos != -2) {
      cout << "build from exist version " << endl;
    }
    outNeighbors.push_back(0);
    myVector<uintE> new_ver(data, data+length);
    edges.push_back(new_ver);
    pos = 0;
    degree = length;
  }

  void append(myVector<uintE> data, uintT version) {
    if (pos == -2) {
      outNeighbors.push_back(0);
      myVector<uintE> empty;
      edges.push_back(empty);
      pos = 0;
      degree = 0;
      // abort();
    }
    if (pos == -1) {
      cout << "append to high vertex " << endl;
      return ;
    }
    if (version <= outNeighbors[outNeighbors.size()-1]) {
      return;
    }
    myVector<uintE> new_ver(data);
    edges.careful_push_back(new_ver);
    outNeighbors.push_back(version);
  }

  inline void empty_step(uintT arg) {
    return ;
  }

  inline void empty_step() {  }

  void prepare() {
    myVector<uintE> data;
    build(data);
  }

  void switch_to(uintT version) {
    if (pos == -1) return;
    if (version >= outNeighbors[pos]) {
      for (; pos<outNeighbors.size()-1; pos++) {
        if (outNeighbors[pos+1]>version) {
          break;
        }
      }
    } else {
      if (pos == 0) return;
      for (; pos > 0; pos--) {
        if (outNeighbors[pos-1] <= version) {
          pos -= 1;
          break;
        }
      }
    }
    degree = edges[pos].size();
  }

  void forward(uintT version) {
    if (pos == -1) return;
    for (; pos<outNeighbors.size()-1; pos++) {
      if (outNeighbors[pos+1]>version) {
        break;
      }
    }
    degree = edges[pos].size();
  }
  void backward(uintT version) {
    if (pos == -1) return;
    for (; pos > 0; pos--) {
      if (outNeighbors[pos-1] <= version) {
        pos -= 1;
        break;
      }
    }
    degree = edges[pos].size();
  }
  void reserve(uintT s) {
    if (pos == -1) {
      outNeighbors.reserve(s);
    }
  }
  void setInNeighbor(uintT j, uintE ngh) {}
  void setOutNeighbor(uintT j, uintE ngh) {}
  uintE getInNeighbor(uintT j) {
    if (pos == -1) {
      return outNeighbors[j];
    } else {
      return edges[pos][j];
      // return 0;
    }
  }
  uintE getOutNeighbor(uintT j) {
    return getInNeighbor(j);
  }

  inline uintT getInDegree() {
    return degree;
  }
  uintT getOutDegree() {
    return getInDegree();
  }

  void print() {
    cout << "print begin " << endl;
    if (is_high_degree_vertex()) {
      for (auto i=0; i<getInDegree(); i++) {
        cout << getInNeighbor(i) << endl;
      }
    } else {
      for (auto i=0; i<outNeighbors.size(); i++) {
        cout << "version :" << outNeighbors[i] ;
        if (pos == i) {
          cout << " <-current version ";
        }
        cout << endl;
        for (auto j=0; j<edges[i].size(); j++) {
          cout << edges[i][j] << endl;
        }
      }
    }
    cout << "print end " << endl;
  }

  bool is_ld_to(uintT ver) {
    if (pos == -1) {
      return false;
    }
    if (pos > 0 && outNeighbors[pos] > ver){
      return true;
    }
    if (pos < outNeighbors.size() - 1 && outNeighbors[pos-1] < ver) {
      return true;
    }
    return false;
  }

  int	find(const uintE &val) {
    if (pos == -1) {
      return outNeighbors.find(val);
    } else {
      return edges[pos].find(val);
    }
  }
  HVertex operator=(const HVertex & other) {
    if (&other != this) {
      outNeighbors = other.outNeighbors;
      pos = other.pos;
      edges = other.edges;
      degree = other.degree;
    }
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<HVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G& g) {
     decode_uncompressed::decodeOutNgh<HVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<HVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<HVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<HVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<HVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<HVertex, F>(this, i, f, bits, tmp1);
  }
};

struct symmetricVertex {
#ifndef WEIGHTED
  myVector<uintE> outNeighbors;
#else
  myVector<intE> outNeighbors;
#endif
  void del() { outNeighbors.clear(); }
#ifndef WEIGHTED
symmetricVertex(myVector<uintE> n)
#else
symmetricVertex(myVector<intE> n)
#endif
: outNeighbors(n) {}

symmetricVertex()
#ifndef WEIGHTED
{ outNeighbors = myVector<uintE>(); }
#else
{ outNeighbors = myVector<intE>(); }
#endif

#ifndef WEIGHTED
  uintE* getInNeighbors () { return outNeighbors.data(); }
  const uintE* getInNeighbors () const { return outNeighbors.data(); }
  uintE* getOutNeighbors () { return outNeighbors.data(); }
  const uintE* getOutNeighbors () const { return outNeighbors.data(); }
  uintE getInNeighbor(uintT j) const { return outNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }

  void setInNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setInNeighbors(uintE* _i) { return ; }
  void setOutNeighbors(uintE* _i) { return ; }
#else
  //weights are stored in the entry after the neighbor ID
  //so size of neighbor list is twice the degree
  intE* getInNeighbors () { return outNeighbors.data(); }
  const intE* getInNeighbors () const { return outNeighbors.data(); }
  intE* getOutNeighbors () { return outNeighbors.data(); }
  const intE* getOutNeighbors () const { return outNeighbors.data(); }
  intE getInNeighbor(intT j) const { return outNeighbors[2*j]; }
  intE getOutNeighbor(intT j) const { return outNeighbors[2*j]; }
  intE getInWeight(intT j) const { return outNeighbors[2*j+1]; }
  intE getOutWeight(intT j) const { return outNeighbors[2*j+1]; }
  void setInNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setInWeight(uintT j, intE wgh) { outNeighbors[2*j+1] = wgh; }
  void setOutWeight(uintT j, intE wgh) { outNeighbors[2*j+1] = wgh; }
  void setInNeighbors(intE* _i) {  }
  void setOutNeighbors(intE* _i) {  }
#endif

  void index_delete(uintT pos) { outNeighbors.index_delete(pos); }
  void index_addtion(uintE data, uintT pos) { outNeighbors.index_addtion(data, pos); }
  void push_back(uintE data) { outNeighbors.push_back(data); }
  void push_back(myVector<uintE>::iterator start, myVector<uintE>::iterator end) { outNeighbors.push_back(start, end); }
  void pop_back() { outNeighbors.pop_back(); }
  bool is_high_degree_vertex() { return true; }
  void switch_to(uintT version) {}
  void backward(uintT version) {}
  void forward(uintT version) {}
  myVector<uintE> getLastVersion() { return outNeighbors; }
  void reserve(uintT s) { outNeighbors.reserve(s); }
  void append(myVector<uintE> data, uintT version) {}
  void print() {
    cout << "print begin" << endl;
    for (auto i=0; i<getInDegree(); i++) {
      cout << getInNeighbor(i) << endl;
    }
    cout << "print end" << endl;
  }

  void prepare() {}

  void empty_step(uintT arg) {
    return ;
  }

  inline void empty_step() { return ; }

  bool is_ld_to(uintT ver) {
    return false;
  }

  uintT getInDegree() const { return outNeighbors.size(); }
  uintT getOutDegree() const { return outNeighbors.size(); }
  void setInDegree(uintT _d) { return;  }
  void setOutDegree(uintT _d) { return;  }
  int	find(const uintE &val) { return outNeighbors.find(val); }
  void flipEdges() {}

  symmetricVertex& operator=(const symmetricVertex & other) {
    
    if (&other != this) {
      outNeighbors = other.outNeighbors;
    }
    return *this;
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<symmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G& g) {
     decode_uncompressed::decodeOutNgh<symmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<symmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<symmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<symmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<symmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<symmetricVertex, F>(this, i, f, bits, tmp1);
  }

};

struct asymmetricVertex {
#ifndef WEIGHTED
  myVector<uintE> inNeighbors, outNeighbors;
#else
  myVector<intE> inNeighbors, outNeighbors;
#endif
  // uintT outDegree;
  // uintT inDegree;
  void del() {inNeighbors.clear(); outNeighbors.clear();}

#ifndef WEIGHTED
asymmetricVertex(myVector<uintE> iN, myVector<uintE> oN)
#else
asymmetricVertex(myVector<intE> iN, myVector<intE> oN)
#endif
: inNeighbors(iN), outNeighbors(oN){}

asymmetricVertex()
#ifndef WEIGHTED
{
	inNeighbors = myVector<uintE>();
	outNeighbors = myVector<uintE>();
}
#else
{	
	inNeighbors = myVector<intE>();
	outNeighbors = myVector<intE>();
}
#endif

#ifndef WEIGHTED
  uintE* getInNeighbors () { return inNeighbors.data(); }
  const uintE* getInNeighbors () const { return inNeighbors.data(); }
  uintE* getOutNeighbors () { return outNeighbors.data(); }
  const uintE* getOutNeighbors () const { return outNeighbors.data(); }
  uintE getInNeighbor(uintT j) const { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  // TODO:: 
  void setInNeighbors(uintE* _i) { return;  }
  void setOutNeighbors(uintE* _i) { return;  }
#else
  intE* getInNeighbors () { return inNeighbors.data(); }
  const intE* getInNeighbors () const { return inNeighbors.data(); }
  intE* getOutNeighbors () { return outNeighbors.data(); }
  const intE* getOutNeighbors () const { return outNeighbors.data(); }
  intE getInNeighbor(uintT j) const { return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) const { return outNeighbors[2*j]; }
  intE getInWeight(uintT j) const { return inNeighbors[2*j+1]; }
  intE getOutWeight(uintT j) const { return outNeighbors[2*j+1]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[2*j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setInWeight(uintT j, uintE wgh) { inNeighbors[2*j+1] = wgh; }
  void setOutWeight(uintT j, uintE wgh) { outNeighbors[2*j+1] = wgh; }
  void setInNeighbors(intE* _i) {  }
  void setOutNeighbors(intE* _i) {  }
#endif

  uintT getInDegree() const { return inNeighbors.size(); }
  uintT getOutDegree() const { return outNeighbors.size(); }
  uintT	find(const uintE &val) { return outNeighbors.find(val); }
  // TODO::::
  void setInDegree(uintT _d) { return;  }
  void setOutDegree(uintT _d) { return;  }
  void flipEdges() { inNeighbors.swap(outNeighbors);  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<asymmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G &g) {
    decode_uncompressed::decodeOutNgh<asymmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_uncompressed::decodeOutNghSparse<asymmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_uncompressed::decodeOutNghSparseSeq<asymmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_uncompressed::copyOutNgh<asymmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_uncompressed::countOutNgh<asymmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_uncompressed::packOutNgh<asymmetricVertex, F>(this, i, f, bits, tmp1);
  }

};

#endif
