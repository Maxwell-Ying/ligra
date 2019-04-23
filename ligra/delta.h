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
#include <get_mem.h>
#include <unordered_map>

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

bool equal(pair<uintE, uintE> a, pair<uintE, uintE> b) {
  return a.first == b.first && a.second == b.second;
}

bool lessthan(pair<uintE, uintE> a, pair<uintE, uintE> b) {
  return a.first < b.first || (a.first == b.first && a.second < b.second);
}

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
  double add_rate = 0.1;
  double delta_rate = 0.001;

  delta_log(myVector<intTriple> _deltalog):deltaLog(_deltalog){}

  delta_log() {}

  delta_log(graph<vertex> &graph, double _add_rate, double _delta_rate, bool diff) {
    add_rate = _add_rate;
    delta_rate = _delta_rate;
    ver = graph.get_version();
    ver_end = ver + 1;
    for(uintT i = 0; i < graph.n; i++) {
      myVector<pair<uintE, intE>> add_and_del = get_delta(graph.getvertex(i), graph.n);
      if (add_and_del.size() > 1 && add_and_del.size()>graph.getvertex(i)->getOutDegree()*delta_rate*2){
        cout << add_and_del.size() << " " << graph.getvertex(i)->getOutDegree() << endl;
      }
      for (auto deledge : add_and_del) {
        deltaLog.push_back(make_pair(i, make_pair(deledge.first, deledge.second)));
      }
    }
    quickSort(deltaLog.data(), deltaLog.size(), logLT());
  }

  myVector<pair<uintE, intE>> get_delta(vertex * v, uintE max) {
    uintT s = v->getInDegree();
    if (s == 0) {
      if (randomFloatBiggerThan(1-100.0/max)) {
        return myVector<pair<uintE, intE>>(1, make_pair(rand()%max, -1));
      }
      return myVector<pair<uintE, intE>>();
    }
    myVector<pair<uintE, intE>> ret;
    double remain = s * delta_rate;

    vector<uintT> add_new, del_new;
    
    while (remain > 0) {
      if (remain < 1 && randomFloatBiggerThan(remain)) break;
      if (randomFloatBiggerThan(1-add_rate)) {
        uintT tmp = rand() % max;
        if (v->find(tmp) == -1) {
          add_new.push_back(tmp);
        } else {
          continue;
        }
      } else {
        del_new.push_back(rand()%s);
      }
      remain -= 1;
    }
    sort(add_new.begin(), add_new.end());
    add_new.erase(unique(add_new.begin(), add_new.end()), add_new.end());
    sort(del_new.begin(), del_new.end());
    del_new.erase(unique(del_new.begin(), del_new.end()), del_new.end());
    for (auto i: add_new) {
      ret.push_back(make_pair(i, -1));
    }
    for (auto i: del_new) {
      ret.push_back(make_pair(v->getInNeighbor(i), i));
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
      intE pos = graph.getvertex(i.first)->find(i.second);
      deltaLog.push_back(make_pair(i.first, make_pair(i.second, pos)));
    }
    quickSort(deltaLog.data(), deltaLog.size(), logLT());
  }

  int size() {
    return deltaLog.size();
  }

  uintT get_next(uintT current) {
    for (uintT i=current+1; i<size(); i++) {
      if (deltaLog[i].first != deltaLog[current].first) {
        return i;
      }
    }
    return size();
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
  
  getline(infile, code);
  delta_log<vertex> d;
  // infile >> version_info;
  getline(infile, version_info);
  cout << "version info " << version_info << endl;
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
      pos = g.getvertex(start)->find(end);
      if (pos == -1) {
        cout << "not found edge " << end << " in vertex " << start << " 's neighbor in version" << d.ver << endl; 
        abort();
      }

      d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
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
  int del_size = -1;
  myVector<uintE> ldVertex;

  delta() {}

  delta(int _vs, int _ve, myVector<uintE>& v, myVector<uintE>& p, myVector<uintE>& dap) :
    vs(_vs), ve(_ve), vertexs(v), positions(p), dstAndPos(dap) {}

  delta(delta_log<vertex> &dlg, graph<vertex> & graph) {
    if (!graph.is_hybrid()) {
      gen(dlg, graph);
    }
    else {
      gen_hybrid(dlg, graph);
    }
    // gen_hybrid(dlg, graph);
    vs = graph.get_version();
    ve = dlg.ver_end;
    del_size = -1;
    del_size = get_del_number();
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
    // cout << "gen finish " << endl;
  }

  void gen_hybrid(delta_log<vertex> &log, graph<vertex> &graph) {
    uintT pos = 0;
    uintT next;
    while (pos < log.size()) {
      next = log.get_next(pos);
      uintE start = log.deltaLog[pos].first;
      vertex * vtx = graph.getvertex(start);
      if (vtx->is_high_degree_vertex()) {
        vertexs.push_back(start);
        positions.push_back(dstAndPos.size());
        bool waitformid = true;
        for (uintT i=pos; i<next; i++) {
          uintE data = log.deltaLog[i].second.first;
          intE weight = log.deltaLog[i].second.second;
          if (weight >= 0) {
            dstAndPos.push_back(data);
            dstAndPos.push_back(weight);
          } else {
            if (waitformid) {
              positions.push_back(dstAndPos.size());
              waitformid = false;
            }
            dstAndPos.push_back(data);
          }
        }
        if (waitformid) {
          positions.push_back(dstAndPos.size());
        }
      } else {
        myVector<uintE> new_ver(vtx->getLastVersion());
        for (auto i=pos; i<next; i++) {
          uintE data = log.deltaLog[i].second.first;
          intE weight = log.deltaLog[i].second.second;
          if (weight >= 0) {
            new_ver.index_delete(weight);
          } else {
            new_ver.push_back(data);
          }
        }
        vtx->append(new_ver, graph.get_version()+1);
        ldVertex.push_back(start);
      }
      pos = next;
    }
    
  }

  delta(graph<vertex> & graph, myVector<delta<vertex> *> path, uintT ver_begin, uintT ver_end) {
    myVector<pair<uintE, uintE>> tmp_add, tmp_del;
    // cout << "into delta gen" << endl;
    for (auto i : path) {
      if (i->vs == ver_begin) {
        i->pack_add_items(&tmp_add);
        i->pack_del_items(&tmp_del);
        ver_begin = i->ve;
      } else if (i->ve == ver_begin) {
        i->pack_add_items(&tmp_del);
        i->pack_del_items(&tmp_add);
        ver_begin = i->vs;
      } else {
        cout << "error in trace path" << endl;
      }
    }

    // 
    quickSort(tmp_add.data(), tmp_add.size(), logLT());
    quickSort(tmp_del.data(), tmp_del.size(), logLT());

    myVector<pair<uintE, uintE>> ret_add, ret_del;
    uintT pa=0, pd=0;
    while (pa < tmp_add.size() && pd < tmp_del.size()) {
      if (equal(tmp_add[pa], tmp_del[pd])) {
        pa ++; pd++;
      } else if (lessthan(tmp_add[pa], tmp_del[pd])) {
        ret_add.push_back(make_pair(tmp_add[pa].first, tmp_add[pa].second));
        pa ++;
      } else {
        ret_del.push_back(make_pair(tmp_del[pd].first, tmp_del[pd].second));
        pd++;
      }
    }
    if (pa == tmp_add.size()) {
      ret_del.push_back(tmp_del.begin() + pd, tmp_del.end());
    } else if (pd == tmp_del.size()) {
      ret_add.push_back(tmp_add.begin() + pa, tmp_add.end());
    }

    delta_log<vertex> dlg = delta_log<vertex>(graph, ret_add, ret_del);
    
    // cout << "before gen " << endl;
    if (!graph.is_hybrid()) {
      gen(dlg, graph);
    }
    else {
      gen_hybrid(dlg, graph);
    }
    // gen_hybrid(dlg, graph);

    vs = graph.get_version();
    ve = get_path_end(path, graph.get_version());
    del_size = -1;
    del_size = get_del_number();
    auto _ld = graph.get_ldvertex(ve);
    ldVertex.push_back(_ld.begin(), _ld.end());
    // cout << ldVertex.size() << endl;
  }

  void write_edge_entry(const char * filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename << endl;
      abort();
    }  
    outfile << "DELTA_FILE" << endl;
    
    outfile << vs << " " << ve << endl;
    outfile << get_add_size() << endl << get_del_number() << endl;

    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT mid, end;
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      for (auto j=mid; j<end; j++) {
        outfile << from << " " << dstAndPos[j] << endl;
      }
    }
    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT start, mid;
      start = positions[2*i];
      mid = positions[2*i+1];
      for (auto j=start; j<mid; j+=2) {
        outfile << from << " " << dstAndPos[j] << endl;
      }
    }
    outfile.close();
  }

  void write_deltalog(const char * filename) {
    ofstream outfile;
    outfile.open(filename);
    if (!outfile.is_open()) {
      cout << "fail to open file " << filename << endl;
      abort();
    }
    outfile << "DELTA_LOG_FILE" << endl;
    outfile << vs << " " << ve << endl;
    outfile << get_add_size() + get_del_number() << endl;
    for (auto i=0; i< vertexs.size(); i++) {
      auto from = vertexs[i];
      uintT start, mid, end;
      start = positions[2*i];
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      for (auto j=start; j<mid; j+=2) {
        outfile << from << " " << dstAndPos[j] << " " << dstAndPos[j+1] << endl;
      }
      for (auto j=mid; j<end; j++) {
        outfile << from << " " << dstAndPos[j] << " " << -1 << endl;
      }
    }
    outfile.close();
  }

  uintT get_path_end(myVector<delta<vertex> *> &path, uintT ver_begin) {
    if (path.empty()) {
      cout << "error in get path end" << endl;
    }
    // cout << ver_begin << " to ";
    for (auto i: path) {
      if (i->vs == ver_begin) {
        ver_begin = i->ve;
      } else if (i->ve == ver_begin) {
        ver_begin = i->vs;
      }
      else {
        cout << "invalid path from " << endl;
        print_path(path);
      }
    }
    // cout << ver_begin << endl;
    return ver_begin;
  }

  // return all add pair of delta

  void pack_add_items(vector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      uintT start, end;
      start = positions[2*i+1];
      if (i == vertexs.size()-1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      if (end > dstAndPos.size()) {
        cout << i << " " << vs << " " << end << " " << start << " " << dstAndPos.size() << endl;
        abort();
      }
      for (auto j=start; j<end; j++) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  void pack_add_items(myVector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      uintT start, end;
      start = positions[2*i+1];
      if (i == vertexs.size()-1) {
        end = dstAndPos.size();
      } else {
        end = positions[2*i+2];
      }
      if (end > dstAndPos.size()) {
        cout << i << " " << vs << " " << end << " " << start << " " << dstAndPos.size() << endl;
        abort();
      }
      for (auto j=start; j<end; j++) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  // return all del pair of delta

  void pack_del_items(vector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      for (auto j=positions[i*2]; j<positions[i*2+1]; j+=2) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  void pack_del_items(myVector<pair<uintE, uintE>> * ret) {
    for (auto i=0; i<vertexs.size(); i++) {
      for (auto j=positions[i*2]; j<positions[i*2+1]; j+=2) {
        ret->push_back(make_pair(vertexs[i], dstAndPos[j]));
      }
    }
  }

  uintE get_delta_size() {
    return dstAndPos.size();
  }

  uintT get_del_number() {
    if (del_size >= 0) {
      return (uintT) del_size;
    }
    uintT ret=0;
    uintT num = positions.size() / 2;

    for (auto i=0; i<num; i++) {
      ret += (positions[2*i+1] - positions[2*i]) / 2;
    }
    del_size = ret;
    return ret;
  }

  uintT get_del_size() {
    return 2* get_del_number();
  }

  uintT get_add_size() {
    return dstAndPos.size() - get_del_size();
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

  int find_edge(uintE s, uintE e) {
    for (auto i=0; i<vertexs.size(); i++) {
      if (vertexs[i] != s) {
        continue;
      }
      uintT begin, mid, end;
      begin = positions[2*i];
      mid = positions[2*i+1];
      if (i == vertexs.size() - 1) {
        end = dstAndPos.size();
      } else {
        end = dstAndPos[2*i+2];
      }
      for (auto j=begin; j<mid; j++) {
        if (dstAndPos[j] == e) {
          return dstAndPos[j+1];
        }
      }
      for (auto j=mid; j<end; j++) {
        if (dstAndPos[j] == e) {
          return -1;
        }
      }
    }
    return -2;
  }

  void print_path(myVector<delta<vertex> *> path) {
    if (path.empty()) {
      cout << "empty path" << endl;
      abort();
      return;
    }
    for (auto i: path) {
      cout << i-> vs << " -> " << i->ve << endl;
    }
    cout << endl;
  }
};

template <class vertex> 
vector<uintT> get_trace(vector<delta<vertex>> &deltas, uintT target) {
  vector<uintT> ret;
  while (target) {
    for (auto i=0; i<deltas.size(); i++) {
      if (deltas[i].ve == target) {
        ret.push_back(i);
        target = deltas[i].vs;
      }
    }
  }
  return ret;
}

template <class vertex> 
vector <uintT> get_revert_path(vector<delta<vertex>> &deltas, uintT from, uintT to) {
  if (from == to) {
    vector<uintT> ret;
    return ret;
  }

  vector<uintT> from_trace = get_trace(deltas, from);
  vector<uintT> to_trace = get_trace(deltas, to);

  while (!from_trace.empty() && !to_trace.empty() 
          && (from_trace[from_trace.size()-1] == (to_trace[to_trace.size()-1]))) {
    from_trace.pop_back();
    to_trace.pop_back();
  }

  reverse(to_trace.begin(), to_trace.end());

  from_trace.insert(from_trace.end(), to_trace.begin(), to_trace.end());
  return from_trace;
}

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
  cout << endl;
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

  vector<delta<vertex> *> get_trace(uintT target) {
    vector<delta<vertex> *> ret;
    bool edit = true;
    while (target && edit) {
      edit = false;
      for (auto i=0; i<edges.size(); i++) {
        if (edges[i].second == target) {
          ret.push_back(get_edge(edges[i].first, edges[i].second));
          target = edges[i].first;
          edit = true;
        }
      }
    }
    return ret;
  }

  myVector <delta<vertex> *> get_revert_path(uintT from, uintT to) {
    if (from == to) {
      myVector<delta<vertex> *> ret;
      return ret;
    }

    vector<delta<vertex> *> from_trace = get_trace(from);
    vector<delta<vertex> *> to_trace = get_trace(to);

    while (!from_trace.empty() && !to_trace.empty() 
          && (from_trace[from_trace.size()-1]->vs == to_trace[to_trace.size()-1]->vs)
          && (from_trace[from_trace.size()-1]->ve == to_trace[to_trace.size()-1]->ve)) {
      from_trace.pop_back();
      to_trace.pop_back();
    }

    reverse(to_trace.begin(), to_trace.end());

    from_trace.insert(from_trace.end(), to_trace.begin(), to_trace.end());
    myVector<delta<vertex> *> tmp;
    for (auto i: from_trace) {
      tmp.push_back(i);
    }
    return tmp;
  }

  // update tree, often just pushback, sometimes add a totally new delta and update chain
  void update_tree(graph<vertex> &graph, delta<vertex> * moved) {
    // try to get a vector contains distance between new version and others. 

    ChronoTimer cter;
    auto current_upper = get_new_upper(moved);
    
    if (current_upper == moved->vs) {
      tree.push_back(moved);
      chain.erase(chain.begin());
    } else {
      auto path = get_revert_path(graph.get_version(), current_upper);
      jump(graph, path);
      
      path = get_revert_path(current_upper, moved->ve);
      
      delta<vertex> *new_delta = new delta<vertex> (graph, path, current_upper, moved->ve);

      edges.push_back(make_pair(current_upper, moved->ve));
      remove_edge(moved->vs, moved->ve);
      cout << "new edge " << current_upper << " => " << moved->ve
           << " instead of " << moved->vs << " -> " << moved->ve << endl;
      delete moved;
      
      tree.push_back(new_delta);
      
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
    if (graph.get_version() != chain[1]->vs) {
      cout << "error in update chain, graph version is " << graph.get_version() << " and chain from " << chain[0]->vs << endl; 
      abort();
    }
    vector<delta<vertex> *> new_chain;

    // chain[0] is previous deleted in update_tree
    for (auto i=1; i<chain.size(); i++) {
      myVector<delta<vertex> *> tmp;
      tmp.push_back(chain[i]);
      auto new_delta = new delta<vertex>(graph, tmp, graph.get_version(), graph.get_version() + 1);

      new_chain.push_back(new_delta);
      delete chain[i];
      apply(graph, *new_delta);
    }
    chain = new_chain;
    // print_edges();
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
    return ret > 0 ? ret : ret * -2;
  }

  // get accurate delta size
  int accurate_delta(myVector<delta<vertex> *> path) {
    vector<pair<uintE, uintE>> ret_add, ret_del;
    for (auto i : path) {
      i->pack_add_items(&ret_add);
      i->pack_del_items(&ret_del);
    }
    sort(ret_add.begin(), ret_add.end(), logLT());
    sort(ret_del.begin(), ret_del.end(), logLT());
    uintT pa = 0, pd = 0;
    uintT count_add=0, count_del=0;
    while (pa < ret_add.size() && pd < ret_del.size()) {
      if (equal(ret_add[pa], ret_del[pd])) {
        pa++; pd++;
      } else if (lessthan(ret_add[pa], ret_del[pd])) {
        count_add++;
        pa++;
      } else {
        pd ++;
        count_del++;
      }
    }
    if (pa < ret_add.size()) {
      count_add += ret_add.size() - pa;
    } else if (pd < ret_del.size()) {
      count_del += ret_del.size() - pd;
    } else if (pa > ret_add.size() || pd > ret_del.size()){
      cout << "something error in accurate_del" << endl;
    }
    
    return count_add + count_del * 2;
  }

  // get where the new edge of graph insert. 
  uintT get_new_upper(delta<vertex> * moved) {
    uintT tree_size = tree.size();

    int current_min = moved->get_delta_size();
    
    uintT current_upper = tree_size;

    // moved.vs == tree_size  no need to loop from tree_size
    for (int i=(int)tree_size-1; i>=0; i--) {
      
      myVector<delta<vertex> *> path = get_revert_path(moved->vs, i);
      if (path.empty()) {
        continue;
      }
      path.push_back(moved);
      auto tmp = estimate_delta(path);
      // cout << "estimate " << moved->vs << " -> " << i << " " << tmp << " " << current_min << endl;
     
      if (tmp < current_min) {
        int new_min = accurate_delta(path);
        // cout << "accurate " << moved->vs << " -> " << i << " " << new_min << " " << current_min << endl;
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
      print_edges();
      abort();
      return NULL;
    }

    for (auto i = 0; i < tree.size(); i++) {
      auto t = tree[i];
      if ((t->vs == b && t->ve == e)) {
        return t;
      }
    }
    for (auto i = 0; i < chain.size(); i++) {
      auto c = chain[i];
      if ((c->vs == b && c->ve == e)) {
        return c;
      }
    }
    cout << "not found edge from " << b << " to " << e << endl;
    print_edges();
    return NULL;
  }

  void append(graph<vertex> & graph, delta<vertex>* d) {
    chain.push_back(d);
    // print_edges();
    edges.push_back(make_pair(d->vs, d->ve));
    repart(graph);
  }

  void repart(graph<vertex> & graph) {
    if (chain.size() > divide) {
      delta<vertex> * moved = chain[0];
      update_tree(graph, moved);
      
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
      cout << i->get_add_size() << " vs " << i->get_del_number() << endl;
    }
    cout << "chain:" << endl;
    for (auto i:chain) {
      cout << i->vs << " -> " << i->ve << endl;
      cout << i->get_add_size() << " vs " << i->get_del_number() << endl;
    }
  }
};

template <class vertex>
struct bigDelta {
  myVector<uintE> versions;
  myVector<uintE> vertexs;
  myVector<uintE> positions;
  myVector<uintE> dstAndPos;

  myVector<uintT> ldVersions;
  myVector<uintE> ldVertex;
  myVector<int> add_number;

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

    ldVersions.push_back(ldVertex.size());
    ldVertex.push_back(da.ldVertex.begin(), da.ldVertex.end());

    int b = (int)da.get_add_size();
    int c = (int)da.get_del_number();
    int a =  b-c ;

    add_number.push_back(a);

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

  int get_add_number(uintT ver) {
    if (ver >= add_number.size()) {
      cout << "ver " << ver << " size " << add_number.size() << endl;
      abort();
    }
    if (add_number.size() != versions.size()) {
      cout << "add " << add_number.size() << " version " << versions.size() << endl;
      abort();
    }
    return add_number[ver];
  }

  int get_add_number(uintT start, uintT end) {
    if (start == end) {
      return 0;
    } else if (start > end) {
      return get_add_number(end, start) * -1;
    }
    
    int ret = 0;
    for (auto i=start; i<end; i++) {
      ret += get_add_number(i);
    }
    
    return ret;
  }
};


template <class vertex>
int apply(graph<vertex> & graph, delta<vertex> & da) {
  // version check
  // print_mem("BellmanFord");
  if (graph.get_version() != da.vs) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    return -1;
  }
  // cout << "apply from " << da.vs << " to " << da.ve << endl;
  
  // for (auto i=0; i<20; i++) 
  //   cout << i << " " << graph.access_vertex(i) << " " << graph.getvertex(i)->getInDegree() << endl;
  if (graph.is_hybrid()) {
    for (auto i : da.ldVertex) {
      // high_degree_vertex dont reponse this call and wait for following delta information
      // graph.getvertex(i)->switch_to(da.ve);
      graph.getvertex(i)->forward(da.ve);
    }
  }
  
  int count = (int)da.vertexs.size();
  // print_mem("BellmanFord");
  for(int i=0; i<count; i++) {
    // get target vertex
    auto vtmp = graph.getvertex(da.vertexs[i]);

    // do the delete first
    for(int j=da.positions[2*i]; j<da.positions[2*i+1]; j+= 2) {
      vtmp->index_delete(da.dstAndPos[j+1]);
    }

    myVector<uintE>::iterator it = da.dstAndPos.begin();
    std::size_t b, e;
      // then do some add
    if (i == count - 1) {
      b = da.positions[2*i+1];
      vtmp->push_back(it + b, da.dstAndPos.end());
    } else {
      b = da.positions[2*i+1];
      e = da.positions[2*i+2];
      vtmp->push_back(it+b, it+e);
    }
  }
  
  graph.set_version(da.ve);
  // graph.update_m();
  int add = ((int)(da.get_add_size()) - (int)(da.get_del_number()));
  graph.add_m(add);
  return 0;
}

template <class vertex>
int revert(graph<vertex> &graph, delta<vertex> &da) {
  // version check
  if (graph.get_version() != da.ve) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    return -1;
  }

  // cout << "revert from " << da.vs << " to " << da.ve << endl;
  // cout << "dirty data before " << graph.accessAllEdges() << endl;
  ChronoTimer cter;
  if (graph.is_hybrid()) {
    for (auto i : da.ldVertex) {
      // high_degree_vertex dont reponse this call and wait for following delta information
      // graph.getvertex(i)->switch_to(da.vs);
      graph.getvertex(i)->backward(da.vs);
    }
  }

  int count = (int) da.vertexs.size();
  // cout << "after ld " << cter.elapsed() << endl;
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
  // cout << "after index add " << cter.elapsed() << endl;
  graph.set_version(da.vs);
  int add = ((int)(da.get_add_size()) - (int)(da.get_del_number()));
  graph.add_m(-1*add);
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

  if (ga.is_hybrid()) {
    uintT b, e;
    b = bda.ldVersions[version_start];
    if (version_end == version_max) {
      e = bda.ldVertex.size();
    } else {
      e = bda.ldVersions[version_end];
    }
    if (e-b > ga.n) {
      for (auto i=b; i<e; i++) {
        // high_degree_vertex dont reponse this call and wait for following delta information
        // ga.getvertex(bda.ldVertex[i])->switch_to(version_end);
        ga.getvertex(bda.ldVertex[i])->forward(version_end);
      }
    } else {
      for (auto i=0; i<ga.n; i++) {
        // ga.getvertex(i)->switch_to(version_end);
        ga.getvertex(i)->forward(version_end);
      }
    }
  }
  
  uintT add_count = 0, del_count = 0;
  for (auto ver = version_start; ver < version_end; ver++) {
    int count = (ver == version_max - 1 ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver]);

    int prefix = bda.versions[ver];
    {parallel_for(int i=prefix; i<prefix + count; i++) {
      vertex* vtmp = ga.getvertex(bda.vertexs[i]);
      
      // do the delete first
      for(int j=bda.positions[2*i]; j<bda.positions[2*i+1]; j+= 2) {
        vtmp->index_delete(bda.dstAndPos[j+1]);
        add_count ++;
      }
      myVector<uintE>::iterator it = bda.dstAndPos.begin();
      std::size_t b, e;
      // then do some add
      if (ver == version_max - 1 && i == prefix + count - 1) {
        b = bda.positions[2*i+1];
        vtmp->push_back(it + b, bda.dstAndPos.end());
        del_count += bda.dstAndPos.size() - b;
      } else {
        b = bda.positions[2*i+1];
        e = bda.positions[2*i+2];
        vtmp->push_back(it+b, it+e);
        del_count += e-b;
      }
    }}
  }
  ga.set_version(version_end);
  // ga.update_m();
  ga.add_m(bda.get_add_number(version_start, version_end));
  return 0;
}

template <class vertex> 
int backward(graph<vertex> &ga, bigDelta<vertex> &bda, int step) {
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
  ChronoTimer cter;
  if (ga.is_hybrid()) {
    uintT b, e;
    b = bda.ldVersions[version_end];
    if (version_start == version_max) {
      e = bda.ldVertex.size();
    } else {
      e = bda.ldVersions[version_start];
    }
    if (e-b < ga.n || step < 20) { 
      for (auto i=b; i<e; i++) {
      // high_degree_vertex dont reponse this call and wait for following delta information
        // ga.getvertex(bda.ldVertex[i])->switch_to(version_end);
        // ga.getvertex(bda.ldVertex[i])->empty_step(0);
        ga.getvertex(bda.ldVertex[i])->backward(version_end);
      }
    } else {
      for (auto i=0; i<ga.n; i++) {
        // ga.getvertex(i)->switch_to(version_end);
        // ga.getvertex(i)->empty_step(0);
        ga.getvertex(i)->backward(version_end);
      }
    }
  }
  // cout << "after ld " << cter.elapsed() << endl;
  
  for (int ver = version_start-1; ver > version_end-1; ver--) {
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
  // cout << "after add and delete" << cter.elapsed() << endl;
  ga.set_version(version_end);

  // ga.update_m();
  int diff_m = bda.get_add_number(version_start, version_end);

  ga.add_m(diff_m);
  // cout << "after add m" << cter.elapsed() << endl;
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, bigDelta<vertex>& bda, int target) {
  if (target < 0 || target > bda.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << bda.get_max_version() << endl;
    return -1;
  }
  cout << "jump " << target << " " << ga.get_version() << endl;
  return forward(ga, bda, target - (int) ga.get_version());
}

template <class vertex>
int jump(graph<vertex>& ga, myVector<delta<vertex> *>&path) {
  // if (!path.empty())
    // print_path<vertex>(path);
  // cout << ga.get_version() << endl;
  // ChronoTimer cter;
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
    // cout << cter.elapsed() << endl;
  }
  return 0;
}

template <class vertex>
int jump(graph<vertex>& ga, versionGraph<vertex> & vg, int target) {
  // first update for initial version.
  // vg.update_vertex_size(ga);
  if (ga.get_version() == target) {
    cout << "jump to same version, nothing happened" << endl;
    return 0;
  } 
  auto path = vg.get_revert_path(ga.get_version(), target);
  cout << "before print path " << endl;
  print_path<vertex>(path);
  cout << "after print path " << endl;
  if (path.empty()) {
    cout << "fail to jump from " << ga.get_version() << " to " << target <<endl;
  }
  jump(ga, path);
  // vg.update_vertex_size(ga);
  return 0;
}

#endif
