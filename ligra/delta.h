#ifndef DELTA_H
#define DELTA_H

#include <random>
#include "graph.h"
#include "myVector.h"
#include "vertex.h"
#include <map>
#include "quickSort.h"
#include <cstdlib>
#include <time.h>
#include "myutil.h"
#include <deque>
#include <set>

typedef pair<uintE, pair<uintE, intE>> intTriple;

// reload of () for sort
struct logLT {
  bool operator () (intTriple a, intTriple b) {
    if (a.first != b.first) {
      return a.first < b.first;
    } else if (a.second.second != b.second.second) {
      return a.second.second > b.second.second;
    } else {
      return a.second.first < b.second.first;
    }
  }
  bool operator () (pair<uintT, uintT> a, pair<uintT, uintT> b) {
    if (a.first != b.first) {
      return a.first < b.first;
    }
    return a.second < b.second;
  }
  // sort by : 1, 2, 3+4, 3
  bool operator () (pair<pair<uintT, uintT>, pair<uintT, uintT>> a, pair<pair<uintT, uintT>, pair<uintT, uintT>> b) {
    if (a.first.first != b.first.first) {
      return a.first.first < b.first.first;
    } else if (a.first.second != b.first.second) {
      return a.first.second < b.first.second;
    } else if (a.second.first + a.second.second != b.second.first + b.second.second) {
      return a.second.first + a.second.first < b.second.first + b.second.second;
    } else {
      return a.second.first < b.second.first;
    }
  }
};



bool randomFloatBiggerThan(double level) {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> ran(0, 1.0);
  return ran(gen) >= level;
}

uintT get_next(myVector<uintT> &arr, uintT target, uintT near) {
  uintT m = arr.size();
  if (near == 0) {
    near = 1;
  }
  while (arr[near] == target) {
    if (near == m-1) {
      return near;
    }
    near ++;
  }
  return near-1;
}

uintT bin_search(myVector<uintT> &arr, uintT target) {
  uintT start = 0;
  uintT end = arr.size() - 1;
  while(end - start > 1) {
    if (arr[end] == target) {
      return get_next(arr, target, end);
    }
    if (arr[start] == target) {
      return get_next(arr, target, start);
    }
    uintT middle = (end+start)/2;
    if (arr[middle] > target) {
      end = middle;
    } else {
      start = middle;
    }
  }
  return start;
}

template <class vertex>
struct delta_log{
  myVector<intTriple> deltaLog;
  uintT ver;
  uintT ver_end;
  double add_rate = 0.9;
  double delta_rate = 0.01;
  // static uintT edgenumber;
  // static vector<uintT> arr;
  delta_log(myVector<intTriple> _deltalog):deltaLog(_deltalog){}

  delta_log() {}

  delta_log(graph<vertex> &graph) {
    ver = graph.get_version();
    ver_end = ver + 1;
    for(uintT i = 0; i < graph.n; i++) {
      myVector<pair<uintE, intE>> add_and_del = get_delta(graph.V[i], graph.n);
      if (add_and_del.size() > 1 && add_and_del.size()>graph.V[i].getOutDegree()*delta_rate*2){
        cout << add_and_del.size() << " " << graph.V[i].getOutDegree() << endl;
      }
      for (auto deledge : add_and_del) {
        deltaLog.push_back(make_pair(i, make_pair(deledge.first, deledge.second)));
      }

    }
    quickSort(deltaLog.data(), deltaLog.size(), logLT());
    cout << deltaLog.size() << endl;
  }

  myVector<pair<uintE, intE>> get_delta(vertex v, uintE max) {
    uintT s = v.getInDegree();
    if (s == 0) {
      if (randomFloatBiggerThan(1-100.0/max)) {
        return myVector<pair<uintE, intE>>(1, make_pair(rand()%max, -1));
      }
      return myVector<pair<uintE, intE>>();
    }
    myVector<pair<uintE, intE>> ret;
    double remain = s * delta_rate;
    intE pos = 0;
    while (remain > 0) {
      if (remain < 1 && randomFloatBiggerThan(remain)) break;
      if (randomFloatBiggerThan(1-add_rate)) {
        uintE e = rand()%max;
        ret.push_back(make_pair(e, -1));
      } else {
        ret.push_back(make_pair(v.getInNeighbor(pos), pos));
        pos += 1;
      }
      remain -= 1;
      // not consider duplicate
    }
    return ret;
  }

  delta_log(graph<vertex> & graph, double _add_rate, double _delta_rate) {
    delta_rate = _delta_rate;
    add_rate = _add_rate;
    // uintT edgesum = graph.get_edge_number();
    uintT edgesum = 54684899;
    uintT number = (uintT) (edgesum * delta_rate);
    myVector<uintT> arr(graph.n+1, 0);
    for(auto i = 1; i <= graph.n; i++) {
      arr[i] = graph.getvertex(i-1)->getInDegree() + arr[i-1];
    }
    cout << arr[100000] << " " << arr[graph.n] << endl;
    uintT i = 0;
    map<uintE, uintT> delcount;
    
    while(i < number) {
      uintE start = bin_search(arr, rand()%graph.n);
      uintE end = rand() % graph.n;

      if (randomFloatBiggerThan(add_rate)) {
        delcount[start] = delcount[start] + 1;
      } else {
        deltaLog.push_back(make_pair(start, make_pair(end, -1)));
      }
      i++;
    }
    
    for(auto const & de : delcount) {
      uintT c = de.second;
      uintT d = graph.getvertex(de.first)->getInDegree();
      if (de.second > d) {
        if (d >= 1) {
          c = 1;
        } else {
          c = 0;
        }
      }
      for (auto j=0; j<c; j++) {
        uintE x = graph.getvertex(de.first)->getInNeighbor(j);
        deltaLog.push_back(make_pair(de.first, make_pair(x, j)));
      }
    }
    quickSort(deltaLog.data(), deltaLog.size(), logLT());
  }

  delta_log(graph<vertex> & graph, myVector<pair<uintE, uintE>> &add, myVector<pair<uintE, uintE>> &del) {
    for (auto i : add) {
      deltaLog.push_back(make_pair(i.first, make_pair(i.second, -1)));
    }
    for (auto i : del) {
      deltaLog.push_back(make_pair(i.first, make_pair(i.second, graph.getvertex(i.first)->find(i.second))));
    }
    quickSort(deltaLog.data(), deltaLog.size(), logLT());
  }

  int size() {
    return deltaLog.size();
  }
};

/*  file structure
 *  1: DELTA_LOG_FILE
 *  2: (base version)
 *  3: (size of deltalog)
 *  4: (first int triple of deltalog)
 *  ......
 *  end of file
 * */
template <class vertex>
int write_deltalog_to_file(delta_log<vertex> d, const char * filename) {
  ofstream outfile;
  outfile.open(filename);
  if (!outfile.is_open()) {
    cout << "fail to open file " << filename<< endl;
    abort();
  }
  outfile << "DELTA_LOG_FILE" << endl;
  outfile << d.ver << " " << d.ver_end << endl;
  outfile << d.deltaLog.size() << endl;
  for (auto const dl : d.deltaLog) {
    outfile << dl.first << " " << dl.second.first << " " << dl.second.second << endl;
  }
  outfile.close();
  return 0;
}

template <class vertex>
delta_log<vertex> load_deltalog_from_file(graph<vertex> & g, const char * filename) {
  ifstream infile;
  infile.open(filename);
  if (!infile.is_open()) {
    cout << "error in open file " << filename << endl;
    abort();
  }
  string code, version_info;
  infile >> code;
  delta_log<vertex> d;
  infile >> version_info;
  auto split = version_info.find(" ");
  if (split != version_info.npos) {
    d.ver_end = atoi(version_info.substr(split).c_str());
    d.ver = atoi(version_info.substr(0, split).c_str());
  } else {
    d.ver = atoi(version_info.c_str());
    d.ver_end = d.ver + 1;
  }
  uintE start, end;
  int pos;
  if (! code.compare("DELTA_LOG_FILE")) {
    uintT length;
    infile >> length;
    for (auto i=0; i < length; i++) {
      infile >> start >> end >> pos;
      d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
    }
  } else if (!code.compare("DELTA_FILE")) {
    myVector<pair<uintE, uintE>> add;
    myVector<pair<uintE, uintE>> del;
    uintT add_length, del_length;
    infile >> add_length;
    infile >> del_length;
    
    for (size_t i=0; i < add_length; i++) {
      infile >> start >> end;
      d.deltaLog.push_back(make_pair(start, make_pair(end, -1)));
    }
    for (size_t i=0; i < del_length; i++) {
      infile >> start >> end;
      pos = g.V[start].find(end);
      if (pos == -1) {
        cout << "not found edge " << end << " in vertex " << start << " 's neighbor in version" << d.ver << endl; 
        abort();
      }
      d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
      // cout << "delete edge from " << start << " to " << end << " at pos " << pos << endl;
    }
    quickSort(d.deltaLog.data(), d.deltaLog.size(), logLT());
    d.add_rate = add_length / (add_length + del_length);
    d.delta_rate = (add_length + del_length) / g.m;
  } else {
    cout << "bad delta file" << filename << " with code" << code << endl;
    abort();
  }
  infile.close();
  return d;
}

template <class vertex>
struct delta {
  int vs, ve; // version start and version end
  myVector<uintE> vertexs;
  myVector<uintE> positions;
  myVector<uintE> dstAndPos;

  delta() {}

  delta(int _vs, int _ve, myVector<uintE>& v, myVector<uintE>& p, myVector<uintE>& dap) :
    vs(_vs), ve(_ve), vertexs(v), positions(p), dstAndPos(dap) {}

  delta(delta_log<vertex> &log, graph<vertex> & graph) {
    gen(log, graph);
  }

  void gen(delta_log<vertex> &log,graph<vertex> &graph) {
    int previous = -1;
    int current;
    int data;
    int weight;
    uintE pos = 0;
    int flag = 0; // -1标志该位置上一次进行了删除条目的添加
                  //  1标志进行了添加条目的添加
    for(int i=0; i<log.size(); i++) {
      current = log.deltaLog[i].first;
      data = log.deltaLog[i].second.first;
      weight = log.deltaLog[i].second.second;

      if (current == previous) {
        if (flag == -1 && weight >= 0) {
          dstAndPos.push_back(data);
          dstAndPos.push_back(weight);
          pos += 2;
        } else if (flag == -1 && weight < 0) {
          positions.push_back(pos);
          flag = 1;
          dstAndPos.push_back(data);
          pos += 1;
        } else if (flag == 1 && weight < 0) {
          dstAndPos.push_back(data);
          pos += 1;
        } else {
          std::cout << "unexcept situation" << std::endl;
        }
      } else {
        vertexs.push_back(current);
        positions.push_back(pos);
        if (flag == -1) {
          positions.push_back(pos);
        }
        if (weight < 0) {
          positions.push_back(pos);
          flag = 1;
          dstAndPos.push_back(data);
          pos += 1;
        } else {
          flag = -1;
          dstAndPos.push_back(data);
          dstAndPos.push_back(weight);
          pos += 2;
        }
      }
      previous = current;
    }
    if (positions.size() % 2 && flag == -1) {
      positions.push_back(pos);
    }
    vs = graph.version;
    ve = vs + 1;
  }

  delta(graph<vertex> & graph, myVector<delta<vertex> *> path) {
    myVector<pair<uintE, uintE>> ret_add, ret_del;
    for (auto i : path) {
      for (auto j : i->get_add_items()) {
        ret_add.push_back(j);
      }
      for (auto j : i->get_del_items()) {
        ret_del.push_back(j);
      }
    }
    delta_log<vertex> dlg = delta_log<vertex>(graph, ret_add, ret_del);
    gen(dlg, graph);
  }

  // return all add pair of delta
  myVector<pair<uintE, uintE>> get_add_items() {
    myVector<pair<uintE, uintE>> ret;
    for (auto i=0; i<vertexs.size(); i++) {
      uintT start, end;
      start = positions[2*i+1];
      if (i == vertexs.size()) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      for (auto j=start; j<end; j++) {
        ret.push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
    return ret;
  }

  // return all del pair of delta
  myVector<pair<uintE, uintE>> get_del_items() {
    myVector<pair<uintE, uintE>> ret;
    for (auto i=0; i<vertexs.size(); i++) {
      for (auto j=positions[i*2]; j<positions[i*2+1]; j+=2) {
        ret.push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
    return ret;
  }

  uintE get_delta_size() {
    return dstAndPos.size();
  }

  uintT get_del_number() {
    uintT ret;
    uintT num = positions.size() / 2;
    for (auto i=0; i<num; i++) {
      ret += (positions[2*i+1] - positions[2*i]) / 2;
    }
    return ret;
  }

  uintT get_del_size() {
    return 2* get_del_number();
  }

  uintT get_add_size() {
    return dstAndPos.size() - get_del_size();
  }

  uintT get_vertex_start(uintT index) {
    return positions[2*index];
  }

  uintT get_vertex_bound(uintT index) {
    if (index == vertexs.size() - 1) {
      return dstAndPos.size();
    } else {
      return positions[2*index + 2];
    }
  }

  uintT get_vertex_wall(uintT index) {
    return positions[index*2+1];
  }

  uintT get_vertex_id(uintT index) {
    return vertexs[index];
  }

  void print_vertex(uintT vertex_id, bool add=false) {
    for (auto i=0; i<vertexs.size(); i++) {
      if (vertexs[i] == vertex_id) {
        uintT s = positions[2*i];
        uintT m = positions[2*i+1];
        uintT e = i==vertexs.size()-1 ? dstAndPos.size() : positions[2*i+2];
        for (auto j=s; j<m; j+=2) {
          cout << dstAndPos[j] << " @ " << dstAndPos[j+1] << endl;
        }
        if (add) {
          for (auto j=m; j<e; j++) {
            cout << dstAndPos[j] << endl;
          }
        }
        break;
      }
    }
  }
};

template<class vertex>
void print_path(myVector<delta<vertex> *> path) {
  if (path.empty()) {
    cout << "empty path" << endl;
    abort();
    return;
  }
  for (auto i: path) {
    cout << i-> vs << " -> " << i->ve << endl;
  }
}

template <class vertex>
struct versionGraph {
  int divide = 5;
  myVector<delta<vertex>*> tree;
  vector<delta<vertex>*> chain;
  vector<pair<uintT, uintT>> edges;

  versionGraph() {}

  versionGraph(float di) {
    if (di > 1) {
      divide = (int) di;
    }
  }

  // update tree, often just pushback, sometimes add a totally new delta and update chain
  void update_tree(graph<vertex> &graph, delta<vertex> * moved) {
    // try to get a vector contains distance between new version and others. 
    auto current_upper = get_new_upper(moved);
    // cout << current_upper << endl;
    if (current_upper == moved->vs) {
      tree.push_back(moved);
    } else {
      delta<vertex> *new_delta = new delta<vertex> (graph, get_path(current_upper, moved->ve));
      edges.push_back(make_pair(current_upper, moved->ve));
      remove_edge(moved->vs, moved->ve);
      delete moved;
      tree.push_back(new_delta);
      auto path = get_path(graph.get_version(), current_upper);
      jump(graph, path);
      apply(graph, *new_delta);
      update_chain(graph);
    }    
  }

  void remove_edge(uintT s, uintT e) {
    for (auto i=edges.begin(); i != edges.end(); i++) {
      if ((i->first == s && i->second == e) || (i->first == e && i->second == s)) {
        edges.erase(i);
        return;
      }
    }
    cout << "fail to erase edge from " << s << " to " << e << endl;
  }

  // update only if new delta can not be simply pushback 
  // loop all the chain and reformat delta 
  void update_chain(graph<vertex> &graph) {
    if (graph.get_version() != chain[0]->vs) {
      cout << "error in update chain, graph version is " << graph.get_version() << " and chain from " << chain[0]->vs << endl; 
      abort();
    }
    vector<delta<vertex> *> new_chain;
    
    for (auto i=0; i<chain.size(); i++) {
      myVector<delta<vertex> *> tmp;
      tmp.push_back(chain[i]);
      auto new_delta = delta<vertex>(graph, tmp);
      new_chain.push_back(&new_delta);
      apply(graph, new_delta);
    }
  }

  // get max version end of all version graph(include tree and chain)
  // version of delta of chain is always front of tree.
  // 获取整个版本图的最大版本。
  uintT get_max_version() {
    vector<uintT> ver_ends;
    for (auto i : chain) {
      ver_ends.push_back(i->ve);
    }
    return *max_element(ver_ends.begin(), ver_ends.end());
  }

  // get neighbors of a version in tree.
  // todo: diff tree and chain
  vector<uintT> get_ngh(uintT ver) {
    vector<uintT> ret;
    for (auto edge: edges) {
      if (edge.first == ver) {
        ret.push_back(edge.second);
      } else if (edge.second == ver) {
        ret.push_back(edge.first);
      }
    }
    return ret;
  }

  // get not accurate delta size
  int estimate_delta(myVector<delta<vertex> *> path) {
    int ret = 0;
    for (auto i: path) {
      ret += i->get_add_size();
      ret -= i->get_del_size();
    }
    return ret > 0 ? ret : ret * 2;
  }

  // get accurate delta size
  int accurate_delta(myVector<delta<vertex> *> path) {
    vector<pair<uintE, uintE>> ret_add, ret_del;
    for (auto i : path) {
      for (auto j : i->get_add_items()) {
        ret_add.push_back(j);
      }
      for (auto j : i->get_del_items()) {
        ret_del.push_back(j);
      }
    }
    sort(ret_add.begin(), ret_add.end());
    sort(ret_del.begin(), ret_del.end());
    auto pa = ret_add.begin();
    auto pd = ret_del.begin();
    while (pa != ret_add.end() && pd != ret_del.end()) {
      if (pa == pd) {
        ret_add.erase(pa);
        ret_del.erase(pd);
      } else if (pa > pd) {
        pd ++;
      } else {
        pa ++;
      }
    }
    return ret_add.size() + ret_del.size() * 2;
  }

  // get where the new edge of graph insert. 
  uintT get_new_upper(delta<vertex> * moved) {
    uintT tree_size = tree.size();
    
    int current_min = moved->get_delta_size();
    uintT current_upper = tree_size;
    // cout << "treesize " << tree_size << " " << current_upper << endl;

    // moved.vs == tree_size  no need to loop from tree_size
    for (int i=(int)tree_size-1; i>=0; i--) {
      myVector<delta<vertex> *> path = get_path(moved->vs, i);
      if (path.empty()) {
        continue;
      }
      path.push_back(moved);
      if (estimate_delta(path) < current_min) {
        int new_min = accurate_delta(path);
        if (new_min < current_min) {
          current_min = new_min;
          current_upper = i;
        }
      }
    }
    return current_upper;
  }

  // get delta by version start and version end; return null if not found
  delta<vertex> * get_edge(uintT b, uintT e) {
    if (b == e) {
      cout << "get error edge from " << b << " to " << e << endl;
      return NULL;
    }

    for (auto i = 0; i < tree.size(); i++) {
      auto t = tree[i];
      if ((t->vs == b && t->ve == e) || t->ve == b && t->vs == e) {
        return t;
      }
    }
    for (auto i = 0; i < chain.size(); i++) {
      auto c = chain[i];
      // cout << c->vs << " " << c-> ve << " " << b << " " << e << endl;
      if ((c->vs == b && c->ve == e) || c->ve == b && c->vs == e) {
        return c;
      }
    }
    cout << "not found edge from " << b << " to " << e << endl;
    return NULL;
  }

  uintT get_halfpath_start(myVector<delta<vertex> *> &path, uintT ver_begin) {
    for (auto i: path) {
      if (i->vs == ver_begin) {
        ver_begin = i->ve;
      } else if (i->ve == ver_begin) {
        ver_begin = i->vs;
      }
    }
    return ver_begin;
  }

  myVector<uintT> get_reached(myVector<delta<vertex> *> &path) {
    myVector<uintT> ret;
    for (auto i : path) {
      auto s = i->vs;
      auto e = i->ve;
      if (find(ret.begin(), ret.end(), s) == ret.end()) {
        ret.push_back(s);
      }
      if (find(ret.begin(), ret.end(), e) == ret.end()) {
        ret.push_back(e);
      }
    }
    return ret;
  }

  void get_path(myVector<delta<vertex> *> &path, uintT ver_begin, uintT ver_end) {
    uintT halfbegin = get_halfpath_start(path, ver_begin);
    // cout << "half begin" << halfbegin << endl;
    auto reached = get_reached(path);
    for (auto i: get_ngh(halfbegin)) {
      if (find(reached.begin(), reached.end(), i) != reached.end()) {
        continue;
      }
      path.push_back(get_edge(halfbegin, i));
      if (i==ver_end) {
        return;
      }
      get_path(path, ver_begin, ver_end);
      if (get_halfpath_start(path, ver_begin) == ver_end) {
        return;
      }
      path.pop_back();
    }
  }

  // get path from ver_begin to ver_end
  myVector<delta<vertex> *> get_path(uintT ver_begin, uintT ver_end) {
    myVector<delta<vertex> *> ret;
    if (ver_begin == ver_end) {
      return ret;
    }
    // cout << "get path from " << ver_begin << " to " << ver_end << endl;

    get_path(ret, ver_begin, ver_end);

    if (ret.empty()) {
      cout << "fail in get path from  " << ver_begin << " to " << ver_end << endl;
    }

    // print_path<vertex>(ret);

    return ret;
  }

  void append(graph<vertex> & graph, delta<vertex>* d) {
    chain.push_back(d);
    // print_edges();
    edges.push_back(make_pair(d->vs, d->ve));
    repart(graph);
  }

  void repart(graph<vertex> & graph) {
    if (chain.size() > divide) {
      // tree.push_back(chain.front());
      delta<vertex> * moved = chain[0];
      update_tree(graph, moved);
      // print_edges();
      chain.erase(chain.begin());
      // print_edges();
      // abort();
    } else {
      cout << "chain not reach limit  " << divide << endl;
    }
  }

  void print_edges() {
    if (edges.empty()) {
      cout << "empty version graph " << endl;
      return;
    }
    cout << "edges:" << endl;
    for (auto i:edges) {
      cout << i.first << " -> " << i.second << endl;
    }
    cout << "tree:" << endl;
    for (auto i:tree) {
      cout << i->vs << " -> " << i->ve << endl;
    }
    cout << "chain:" << endl;
    for (auto i:chain) {
      cout << i->vs << " -> " << i->ve << endl;
    }
  }
};

template <class vertex>
struct bigDelta {
  myVector<uintE> versions;
  myVector<uintE> vertexs;
  myVector<uintE> positions;
  myVector<uintE> dstAndPos;
  int append(delta<vertex> da) {
    versions.push_back(vertexs.size());
    vertexs.reserve(vertexs.size()+da.vertexs.size());
    
    for (auto i = 0; i < da.vertexs.size(); i++) {
      vertexs.push_back(da.vertexs[i]);
    }
    
    positions.reserve(positions.size() + da.positions.size());
    int last_length = dstAndPos.size();
    // cout << positions.check_increase() << endl;
    for (auto i = 0; i < da.positions.size(); i++) {
      positions.push_back(da.positions[i]+last_length);
    }
    // cout << positions.check_increase() << endl;
    dstAndPos.reserve(last_length + da.dstAndPos.size());
    for (auto i = 0; i < da.dstAndPos.size(); i++) {
      dstAndPos.push_back(da.dstAndPos[i]);
    }

    return 0;
  }
  int get_max_version() {
    return versions.size();
  }
  uintT size() {
    return dstAndPos.size();
  }
  uintT capacity() {
    return dstAndPos.capacity();
  }
};

template <class vertex>
int apply(graph<vertex> & graph, delta<vertex> & da) {
  // version check
  if (graph.version != da.vs) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    return -1;
  }

  int count = (int)da.vertexs.size();
  
  for(int i=0; i<count; i++) {
    // get target vertex
    auto vtmp = graph.getvertex(da.vertexs[i]);
    // do the delete first
    // cout << da.vertexs[i] << endl;
    for(int j=da.positions[2*i]; j<da.positions[2*i+1]; j+= 2) {
      vtmp->index_delete(da.dstAndPos[j+1]);
    }
    if (i == count - 1) {
      for(int j=da.positions[2*i+1]; j < da.dstAndPos.size(); j+= 1) {
        vtmp->push_back(da.dstAndPos[j]);
      }
    }else {
      for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
        vtmp->push_back(da.dstAndPos[j]);
      }
    }
  }
  graph.version = da.ve;
  graph.update_m();
  return 0;
}

template <class vertex>
int revert(graph<vertex> &graph, delta<vertex> &da) {
  // version check
  if (graph.version != da.ve) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    return -1;
  }

  int count = (int) da.vertexs.size();
  // cout << count << " vertexs" << endl;
  // cout << da.positions.size() << " positions" << endl;
  // cout << da.dstAndPos.size() << " > " << da.positions[da.positions.size()-1] << endl;
  // cout << graph.n << " vertexs" << endl;
  {parallel_for(int i=0; i<count; i++) {
    // {for(int i=0; i<count; i++) {
    
    vertex* vtmp = graph.getvertex(da.vertexs[i]);

    if (i == count - 1) {
      for(int j=da.positions[2*i+1]; j<da.dstAndPos.size(); j+= 1) {
        vtmp->pop_back();
      }
    } else {
      for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
        vtmp->pop_back();
      }
    }
    
    for(int j=da.positions[2*i+1]-2; j>=(int) da.positions[2*i]; j-= 2) {
      vtmp->index_addtion(da.dstAndPos[j], da.dstAndPos[j+1]);
    }
  }}
  // cout << "revert finish" << endl;
  graph.version = da.vs;
  graph.update_m();
  return 0;
}

template <class vertex> 
int forward(graph<vertex> &ga, bigDelta<vertex> &bda, int step = 1) {
  if (step == 0) {
    cout << "forward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return backward(ga, bda, -1 * step);
  }

  int version_start = ga.get_version();
  int version_end = version_start + step;
  int version_max = bda.get_max_version();
  if (version_end > version_max + 1) {
    cout << "try to get far version than exist, nothing happened." << endl;
    cout << version_end << " > " << bda.get_max_version() << endl;
    return -1;
  }
  // cout << "version " << version_start << " " << version_max << endl;
  // cout << ga.V[2165667].outNeighbors.size() << endl;
  
  for (auto ver = version_start; ver < version_end; ver++) {
    int count = (ver == version_max - 1 ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver]);

    int prefix = bda.versions[ver];
    // cout << count << " " << prefix << endl;
    {parallel_for(int i=prefix; i<prefix + count; i++) {
      if (bda.vertexs[i] >= ga.n) {
        ga.V.resize(bda.vertexs[i] + 1);
        ga.n = ga.V.size();
      }
      vertex* vtmp = ga.getvertex(bda.vertexs[i]);
      
      // do the delete first
      for(int j=bda.positions[2*i]; j<bda.positions[2*i+1]; j+= 2) {
        vtmp->index_delete(bda.dstAndPos[j+1]);
      }
      myVector<uintE>::iterator it = bda.dstAndPos.begin();
      std::size_t b, e;
      // then do some add
      if (ver == version_max - 1 && i == prefix + count - 1) {
        // for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
        //   vtmp.outNeighbors.push_back(bda.dstAndPos[j]);
        // }
        b = bda.positions[2*i+1];
        vtmp->push_back(it + b, bda.dstAndPos.end());
      } else {
        // for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
        //   vtmp.outNeighbors.push_back(bda.dstAndPos[j]);
        // }
        b = bda.positions[2*i+1];
        e = bda.positions[2*i+2];
        vtmp->push_back(it+b, it+e);
      }
    }}
  }
  
  ga.version = version_end;
  ga.update_m();
  return 0;
}

template <class vertex> 
int backward(graph<vertex> &ga, bigDelta<vertex> &bda, int step = 1) {
  if (step == 0) {
    cout << "backward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return forward(ga, bda, -1 * step);
  }
  int version_start = ga.get_version();
  if (version_start > bda.get_max_version() + 1) {
    cout << "try to travel from some version not exists: " << bda.get_max_version() << endl;
    return -1;
  }
  int version_end = version_start - step;
  if (version_end < 0) {
    cout << "try to get negative version " << version_end << ", nothing happened." << endl;
    return -1;
  }
  int version_max = bda.get_max_version();

  for (auto ver = version_start-1; ver > version_end-1; ver--) {
    int count = (ver == version_max - 1) ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver];
    
    int prefix = bda.versions[ver];

    {parallel_for(int i=prefix; i<count+prefix; i++) {
      vertex* vtmp = ga.getvertex(bda.vertexs[i]);
      if (ver == version_max - 1 && i == prefix+count-1) {
        for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
          vtmp->pop_back();
        }
      }else {
        for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
          vtmp->pop_back();
        }
      }
      for(int j=bda.positions[2*i+1]-2; j >= (int) bda.positions[2*i]; j-= 2) {
        vtmp->index_addtion(bda.dstAndPos[j], bda.dstAndPos[j+1]);
      }
    }}
  }
  ga.version = version_end;
  ga.update_m();
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, bigDelta<vertex>& bda, int target) {
  if (target < 0 || target > bda.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << bda.get_max_version() << endl;
    return -1;
  }
  return forward(ga, bda, target - ga.get_version());
}

template <class vertex>
int jump(graph<vertex>& ga, myVector<delta<vertex> *>&path) {
  if (!path.empty())
    // print_path<vertex>(path);
  // cout << ga.get_version() << endl;
  for (auto i=0; i<path.size(); i++) {
    auto d = path[i];
    // cout << d->vs << " to " << d->ve << " with graph " << ga.get_version() << endl;
    if (d->vs == ga.get_version()) {
      apply(ga, *d);
    } else if (d->ve == ga.get_version()) {
      revert(ga, *d);
    } else {
      cout << "error in jump with path" << endl;
    }
  }
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, versionGraph<vertex> & vg, int target) {
  // first update for initial version.
  // vg.update_vertex_size(ga);
  auto path = vg.get_path(ga.get_version(), target);
  // print_path<vertex>(path);
  jump(ga, path);
  // vg.update_vertex_size(ga);
  return 0;
}

#endif
