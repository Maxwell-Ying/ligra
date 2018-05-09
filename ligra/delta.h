#ifndef DELTA_H
#define DELTA_H

#include "graph.h"
#include "myVector.h"
#include "vertex.h"

#include "quickSort.h"

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
struct delta {
  int vs, ve; // version start and version end
  myVector<uintE> vertexs;
  myVector<int> positions;
  myVector<uintE> dstAndPos;

  delta(int vs, int ve, myVector<uintE>& v, myVector<int>& p, myVector<uintE>& dap) :
    vs(vs), ve(ve), vertexs(v), positions(p), dstAndPos(dap) {}

  delta(graph<vertex> &graph, myVector<intTriple> log) {
    int add_count =0, del_count = 0;
    int count = log.size();
    for(int i=0; i < count; i++) {
      intTriple tmp = log[i];
      if (tmp.second.second == -1) {
        del_count += 1;
        int pos = graph.V[tmp.first].find(tmp.second.first);
        // 为了编译之便，需要之后的更改。
        
        if (pos == -1) {
          std::cout << "didnt find delete vertex" << std::endl;
        }
        log[i].second.second = pos;
      } else if (tmp.second.second == 1){
        add_count += 1;
        log[i].second.second = -2;
      } else {
        std::cout << "strange weight : " << tmp.second.second << endl;
      }
    }
    // std::cout << log[1].second.second << std::endl;
               // TODO :: 缺乏必要的错误处理，例如在检测到错误条目的时候删除对应的log
                //         假设已经进行处理，之后的数据为正常数据。
    quickSort(log.data(), add_count+del_count, logLT());
    // sort(log);  //  TODO:: 排序尝试使用quicksort中的排序方法，
                //具体用法有待研究
    // std::cout << log[1].second.second << std::endl;

    int previous = -1;
    int current;
    int data;
    int weight;
    int pos = 0;
    int flag = 0; // -1标志该位置上一次进行了删除条目的添加
                  //  1标志进行了添加条目的添加
    for(int i=0; i<add_count+del_count; i++) {
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
  graph.version = da.ve;
  return 0;
}

#endif
