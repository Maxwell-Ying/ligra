#! /usr/bin/env bash

project_home=/public/home/tangwei/project/ligra

graph_data=/public/home/tangwei/graphdata/csr.livejournal-links

result_path=${project_home}/result

factors=(0.001 0.002 0.005 0.01 0.02 0.05 0.1)
times=(200 100 50 50 50 50 30)

for ((i=0;i<7;i++)) ;
do
  result_file=${result_path}/vgraph-f${factors[$i]}-n100.txt
  LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${project_home}/bin/PageRank -f ${factors[$i]} -n 100 -s ${graph_data} > ${result_file}
  sleep 10
done
