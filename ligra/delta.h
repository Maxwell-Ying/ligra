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
typedef pair<uintE, pair<uintE, intE>> intTriple;

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
  double add_rate = 0.9;
  double delta_rate = 0.01;
  static uintT edgenumber;
  static vector<uintT> arr;
  delta_log(myVector<intTriple> _deltalog):deltaLog(_deltalog){}

  delta_log() {}

  delta_log(graph<vertex> &graph) {
    ver = graph.getversion();
    for(uintT i = 0; i < graph.n; i++)
    {
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
    uintT s = v.outNeighbors.size();
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
        ret.push_back(make_pair(v.outNeighbors[pos], pos));
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
    for(auto i = 1; i <= graph.n; i++)
    {
      arr[i] = graph.V[i-1].outNeighbors.size() + arr[i-1];
    }
    cout << arr[100000] << " " << arr[graph.n] << endl;
    uintT i = 0;
    map<uintE, uintT> delcount;
    
    while(i < number){
      uintE start = bin_search(arr, rand()%graph.n);
      uintE end = rand() % graph.n;

      if (randomFloatBiggerThan(add_rate)) {
        delcount[start] = delcount[start] + 1;
      } else {
        deltaLog.push_back(make_pair(start, make_pair(end, -1)));
      }
      i++;
    }
    
    for(auto const & de : delcount)
    {
      uintT c = de.second;
      if (de.second > graph.V[de.first].outNeighbors.size()) {
        if (graph.V[de.first].outNeighbors.size() >= 1) {
          c = 1;
        } else {
          c = 0;
        }
      }
      for (auto j=0; j<c; j++) {
        uintE x = graph.V[de.first].outNeighbors[j];
        deltaLog.push_back(make_pair(de.first, make_pair(x, j)));
      }
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
  outfile << d.ver << endl;
  outfile << d.deltaLog.size() << endl;
  for (auto const dl : d.deltaLog) {
    outfile << dl.first << " " << dl.second.first << " " << dl.second.second << endl;
  }
  outfile.close();
  return 0;
}

template <class vertex>
delta_log<vertex> load_deltalog_from_file(const char * filename) {
  ifstream infile;
  infile.open(filename);
  if (!infile.is_open()) {
    cout << "error in open file " << filename << endl;
    abort();
  }
  string code;
  infile >> code;
  if (code.compare("DELTA_LOG_FILE")) {
    cout << "bad delta file" << filename << " with code" << code << endl;
    infile.close();
    abort();
  }
  delta_log<vertex> d;
  infile >> d.ver;
  uintT length;
  uintT start, end;
  intT pos;
  infile >> length;
  for (auto i=0; i < length; i++) {
    infile >> start >> end >> pos;
    d.deltaLog.push_back(make_pair(start, make_pair(end, pos)));
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

  delta(int _vs, int _ve, myVector<uintE>& v, myVector<uintE>& p, myVector<uintE>& dap) :
    vs(_vs), ve(_ve), vertexs(v), positions(p), dstAndPos(dap) {}

  delta(delta_log<vertex> &log,graph<vertex> &graph) {
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

};

template <class vertex> 
struct deltaVector {
  myVector<delta<vertex>> allDelta;
  // version start from 0, max version should be count + 1
  int get_max_version() {
    return allDelta.size();
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
    vertex& vtmp = graph.getvertex()[da.vertexs[i]];
    // do the delete first
    // cout << da.vertexs[i] << endl;
    for(int j=da.positions[2*i]; j<da.positions[2*i+1]; j+= 2) {
      vtmp.outNeighbors.index_delete(da.dstAndPos[j+1]);
    }
    if (i == count - 1) {
      for(int j=da.positions[2*i+1]; j < da.dstAndPos.size(); j+= 1) {
        vtmp.outNeighbors.push_back(da.dstAndPos[j]);
      }
    }else {
      for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
        vtmp.outNeighbors.push_back(da.dstAndPos[j]);
      }
    }
  }
  graph.version = da.ve;
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
  {parallel_for(int i=0; i<count; i++) {
    vertex& vtmp = graph.getvertex()[da.vertexs[i]];
    if (i == count - 1) {
      for(int j=da.positions[2*i+1]; da.dstAndPos.size(); j+= 1) {
        vtmp.outNeighbors.pop_back();
      }
    } else {
      for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
        vtmp.outNeighbors.pop_back();
      }
    }
    for(int j=da.positions[2*i+1]-2; j>=da.positions[2*i]; j-= 2) {
      vtmp.outNeighbors.index_addtion(da.dstAndPos[j], da.dstAndPos[j+1]);
    }
  }}
  graph.version = da.vs;
  return 0;
}

template <class vertex> 
int forward(graph<vertex> & ga, deltaVector<vertex> & das, int step = 1) {
  if (step == 0) {
    cout << "forward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return backward(ga, das, -1 * step);
  }
  int version_start = ga.getversion();
  int version_end = version_start + step;
  if (version_end > das.get_max_version()+1) {
    cout << "try to get far version than exist, nothing happened." << endl;
    cout << version_end << " > " << das.get_max_version() << endl;
    return -1;
  }
  for (auto i = version_start; i < version_end; i++) {
    apply(ga,das.allDelta[i]);
  }
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
  
  int version_start = ga.getversion();
  int version_end = version_start + step;
  int version_max = bda.get_max_version();
  if (version_end > version_max + 1) {
    cout << "try to get far version than exist, nothing happened." << endl;
    cout << version_end << " > " << bda.get_max_version() << endl;
    return -1;
  }
  for (auto ver = version_start; ver < version_end; ver++) {
    int count = ver == version_max - 1 ? 
                  bda.vertexs.size() - bda.versions[ver] :
                  bda.versions[ver+1] - bda.versions[ver];
    
    int prefix = bda.versions[ver];
    {parallel_for(int i=prefix; i<prefix + count; i++) {
      vertex& vtmp = ga.getvertex()[bda.vertexs[i]];
      // do the delete first
      for(int j=bda.positions[2*i]; j<bda.positions[2*i+1]; j+= 2) {
        vtmp.outNeighbors.index_delete(bda.dstAndPos[j+1]);
      }
      // then do some add
      if (ver == version_max - 1 && i == prefix + count - 1) {
        for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
          vtmp.outNeighbors.push_back(bda.dstAndPos[j]);
        }
      } else {
        for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
          vtmp.outNeighbors.push_back(bda.dstAndPos[j]);
        }
      }
    }}
  }
  ga.version = version_end;
  return 0;
}

template <class vertex> 
int backward(graph<vertex> & ga, deltaVector<vertex> & das, int step = 1) {
  if (step == 0) {
    cout << "backward 0 step, that's to say, nothing happened" << endl;
    return 0;
  }else if (step < 0) {
    return forward(ga, das, -1 * step);
  }
  int version_start = ga.getversion();
  if (version_start > das.get_max_version()) {
    cout << "try to travel from some version not exists: " << das.get_max_version() << endl;
    return -1;
  }
  int version_end = version_start - step;
  if (version_end < 0) {
    cout << "try to get negative version " << version_end << ", nothing happened." << endl;
    return -1;
  }

  for (auto i = version_start; i > version_end; i--) {
    revert(ga,das.allDelta[i-1]);
  }
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
  int version_start = ga.getversion();
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
      vertex& vtmp = ga.getvertex()[bda.vertexs[i]];
      if (ver == version_max - 1 && i == prefix+count-1) {
        for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
          vtmp.outNeighbors.pop_back();
        }
      }else {
        for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
          vtmp.outNeighbors.pop_back();
        }
      }
      for(int j=bda.positions[2*i+1]-2; j >= (int) bda.positions[2*i]; j-= 2) {
        vtmp.outNeighbors.index_addtion(bda.dstAndPos[j], bda.dstAndPos[j+1]);
      }
    }}
  }
  ga.version = version_end;
  return 0;
}

template <class vertex> 
int jump(graph<vertex> & ga, deltaVector<vertex> & das, int target) {
  if (target < 0 || target > das.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << das.get_max_version() << endl;
    return -1;
  }
  return forward(ga, das, target-ga.getversion());
}

template <class vertex>
int jump(graph<vertex>& ga, bigDelta<vertex>& bda, int target) {
  if (target < 0 || target > bda.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << bda.get_max_version() << endl;
    return -1;
  }
  return forward(ga, bda, target - ga.getversion());
}

#endif
