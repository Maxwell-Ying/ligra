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
  };
};

template <class vertex>
struct delta_log{
  myVector<intTriple> deltaLog;
  delta_log(myVector<intTriple> _deltalog):deltaLog(_deltalog){}
  delta_log(graph<vertex> & graph){
    double delta_rate = 0.001;
    uintE delta_number = graph.m * delta_rate;
    double add_rate = 0.6;
    double delete_rate = 1 - add_rate;
    uintE add_number = delta_number * add_rate;
    uintE delete_number = delta_number * delete_rate;
    uintE n = graph.n;
    uintE i = 0;
    uintE j = 0;
    std::map<uintE, myVector<uintE>> log_map;//用来存已在deltalog的边，去重
    while(1){
      uintE vertex_start = rand() % n;
      uintE vertex_end = rand() % n;
      std::map<uintE, myVector<uintE>>::iterator key = log_map.find(vertex_start);
      if(key != log_map.end() && key->second.find(vertex_end) != -1)
        continue;
      else if(key == log_map.end()){
        myVector<uintE> mv;
        mv.push_back(vertex_end);
        log_map[vertex_start] = mv;
        intE pos = graph.V[vertex_start].find[vertex_end];
        intTriple _log;
        _log.first = vertex_start;
        _log.second.first = vertex_end;
        _log.second.second = pos;
        deltaLog.push_back(_log);
        if(pos == -1)
          i++;     
        else
          j++; 
      }
      else if(key != log_map.end() && key->second.find(vertex_end) == -1){
        key->second .push_back(vertex_end);
        intE pos = graph.V[vertex_start].find[vertex_end];
        intTriple _log;
        _log.first = vertex_start;
        _log.second.first = vertex_end;
        _log.second.second = pos;
        deltaLog.push_back(_log);
        if(pos == -1)
          i++;     
        else
          j++;
      }
      if(i == add_number)//默认加边先满
        break;  
    }
    //在图中找减边存到deltalog中
    while(1){
      uintE vertex_start = rand() % n;
      for(i = 0; i < graph.V[vertex_start].getOutDegree; i++){
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
          j++;
        }
        else if(key != log_map.end() && key->second.find(vertex_end) == -1){
          key->second .push_back(vertex_end);
          intTriple _log;
          _log.first = vertex_start;
          _log.second.first = vertex_end;
          _log.second.second = i;
          deltaLog.push_back(_log);
          j++;
        }
      }
      if(j == delete_number)
          break;
    }
    quickSort(deltaLog.data(), add_number + delete_number, logLT());
  }
};
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
template <class vertex>
struct delta {
  int vs, ve; // version start and version end
  myVector<uintE> vertexs;
  myVector<int> positions;
  myVector<uintE> dstAndPos;

  delta(int vs, int ve, myVector<uintE>& v, myVector<int>& p, myVector<uintE>& dap) :
    vs(vs), ve(ve), vertexs(v), positions(p), dstAndPos(dap) {}

  delta(delta_log<vertex> &log,graph<vertex> &graph) {

    int previous = -1;
    int current;
    int data;
    int weight;
    int pos = 0;
    int flag = 0; // -1标志该位置上一次进行了删除条目的添加
                  //  1标志进行了添加条目的添加
    for(int i=0; i<log.size(); i++) {
      current = log[i].first;
      data = log[i].second.first;
      weight = log[i].second.second;

      if (current == previous) {
        if (flag == -1 && weight >= 0) {
          dstAndPos.push_back(data);
          dstAndPos.push_back(weight);
          pos += 2;
        } else if (flag == -1 && weight == -2) {
          positions.push_back(pos);
          flag = 1;
          dstAndPos.push_back(data);
          pos += 1;
        } else if (flag == 1 && weight == -2) {
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
    positions.push_back(pos);
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
int apply(graph<vertex> & graph, delta<vertex> & da) {
  // version check
  if (graph.version != da.vs) {
    std::cout << "version mismatch, delta apply failed, nothing changed" << endl;
    return -1;
  }

  int count = (int)da.vertexs.size();
  // cout << count << endl;
  
  parallel_for(int i=0; i<count; i++) {
    // get target vertex
    vertex& vtmp = graph.getvertex()[da.vertexs[i]];
    // do the delete first
    for(int j=da.positions[2*i]; j<da.positions[2*i+1]; j+= 2) {
      vtmp.outNeighbors.index_delete(da.dstAndPos[j+1]);
    }
    // then do some add
    for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
      vtmp.outNeighbors.push_back(da.dstAndPos[j]);
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
  parallel_for(int i=0; i<count; i++) {
    vertex& vtmp = graph.getvertex()[da.vertexs[i]];
    for(int j=da.positions[2*i+1]; j<da.positions[2*i+2]; j+= 1) {
      vtmp.outNeighbors.pop_back();
    }
    for(int j=da.positions[2*i+1]-2; j>=da.positions[2*i]; j-= 2) {
      vtmp.outNeighbors.index_addtion(da.dstAndPos[j], da.dstAndPos[j+1]);
    }
  }
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
  if (version_end > das.get_max_version()) {
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

#endif
