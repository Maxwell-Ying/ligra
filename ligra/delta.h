#ifndef DELTA_H
#define DELTA_H

#include "graph.h"
#include "myVector.h"
#include "vertex.h"
#include <map>
#include "quickSort.h"
#include <cstdlib>
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

template <class vertex>
struct delta_log{
  myVector<intTriple> deltaLog;
  delta_log(myVector<intTriple> _deltalog):deltaLog(_deltalog){}
  delta_log(graph<vertex> & graph){
    double delta_rate = 0.001;
    uintE delta_number = graph.m * delta_rate;
    double add_rate = 0.6;
    uintE add_number = delta_number * add_rate;
    uintE delete_number = delta_number - add_number;
    uintE n = graph.n;
    uintE add_count = 0;
    uintE del_count = 0;
    std::map<uintE, myVector<uintE>> log_map;//用来存已在deltalog的边，去重
    while(add_count < add_number){  //默认加边先满
      if(del_count >= delete_number) {
        break;
      }
      uintE vertex_start = rand() % n;
      uintE vertex_end = rand() % n;
      std::map<uintE, myVector<uintE>>::iterator key = log_map.find(vertex_start);
      if(key != log_map.end() && key->second.find(vertex_end) != -1)
        continue;
      else if(key == log_map.end()){
        myVector<uintE> mv;
        mv.push_back(vertex_end);
        log_map[vertex_start] = mv;
        intE pos = graph.V[vertex_start].find(vertex_end);
        intTriple _log;
        _log.first = vertex_start;
        _log.second.first = vertex_end;
        _log.second.second = pos;
        deltaLog.push_back(_log);
        if(pos == -1)
          add_count++;     
        else
          del_count++; 
      }
      else if(key != log_map.end() && key->second.find(vertex_end) == -1){
        key->second.push_back(vertex_end);
        intE pos = graph.V[vertex_start].find(vertex_end);
        intTriple _log;
        _log.first = vertex_start;
        _log.second.first = vertex_end;
        _log.second.second = pos;
        deltaLog.push_back(_log);
        if(pos == -1)
          add_count++;     
        else
          del_count++;
      }
    }
    //在图中找减边存到deltalog中
    while(del_count < delete_number){
      uintE vertex_start = rand() % n;
      for(int i = 0; i < graph.V[vertex_start].getOutDegree(); i++){
        uintE vertex_end = graph.V[vertex_start].outNeighbors[i];
        std::map<uintE, myVector<uintE>>::iterator key = log_map.find(vertex_start);
        if(key != log_map.end() && key->second.find(vertex_end) != -1)
          continue;
        else if(key == log_map.end()){
          myVector<uintE> mv;
          mv.push_back(vertex_end);
          log_map[vertex_start] = mv;
          intTriple _log;
          _log.first = vertex_start;
          _log.second.first = vertex_end;
          _log.second.second = i;
          deltaLog.push_back(_log);
          del_count++;
        }
        else if(key != log_map.end() && key->second.find(vertex_end) == -1){
          key->second .push_back(vertex_end);
          intTriple _log;
          _log.first = vertex_start;
          _log.second.first = vertex_end;
          _log.second.second = i;
          deltaLog.push_back(_log);
          del_count++;
        }
        if (del_count >= delete_number) {
          break;
        }
      }
    }
    
    quickSort(deltaLog.data(), add_number + delete_number, logLT());
    
  }
  int size() {
    return deltaLog.size();
  }
};


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
    int count = ver == version_max ? 
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
      if (ver == version_max && i == prefix + count - 1) {
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
    if (ver == version_max && i == prefix+count-1) {
      for(int j=bda.positions[2*i+1]; j<bda.dstAndPos.size(); j+= 1) {
        vtmp.outNeighbors.pop_back();
      }
    }else {
      for(int j=bda.positions[2*i+1]; j<bda.positions[2*i+2]; j+= 1) {
        vtmp.outNeighbors.pop_back();
      }
    }
    
    for(int j=bda.positions[2*i+1]-2; j>=bda.positions[2*i]; j-= 2) {
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
  return forward(ga, das, target-das.get_max_version());
}

template <class vertex>
int jump(graph<vertex>& ga, bigDelta<vertex>& bda, int target) {
  if (target < 0 || target > bda.get_max_version()) {
    cout << "try to jump to version not exist, nothing happened" << endl;
    cout << target << "   " << bda.get_max_version() << endl;
    return -1;
  }
  return forward(ga, bda, target - bda.get_max_version());
}

#endif
