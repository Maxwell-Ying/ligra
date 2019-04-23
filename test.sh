#!/usr/bin/env bash

datasets=(dewiki livejournal enwiki orkut flickr twitter)

dir_prefix=(link-dynamic-dewiki livejournal enwiki orkut flickr twitter)
csr_suffix=(static-dewiki livejournal-links edit-enwiki orkut-links flickr-growth twitter_rv.net)
delta_type=(small normal smalltree delete)

day=$(date +"%m-%d")
result_dir=/public/home/tangwei/project/ligra/result/
target_result_path=${result_dir}/test-${day}
final_result_path=$target_result_path
# if [ -d ${target_result_path} ]
# then
#   for ((j=0; j<5; j++))
#   do
#     final_result_path=${target_result_path}-${j}
#     if [ -d ${final_result_path} ]
#     then
#       sleep 1s
#     else
#       mkdir -p $final_result_path 
#       break
#     fi
#   done
# else 
#   mkdir -p $final_result_path
# fi

delta_num=65

exe=/public/home/tangwei/project/ligra/bin/BellmanFord


# for i in 0 1 3 4
for i in 2 5
do 
  for j in 0 
  do
    data=/public/home/tangwei/graphdata/${dir_prefix[$i]}/csr.${csr_suffix[$i]}
    LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -type ${j} -s ${data} > ${final_result_path}/${datasets[$i]}_basic_${delta_type[j]}.txt
    # LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -type ${j} -g -s ${data} > ${final_result_path}/${datasets[$i]}_tree_${delta_type[j]}.txt
    # LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -type ${j} -h -s ${data} > ${final_result_path}/${datasets[$i]}_hybrid_${delta_type[j]}.txt
    # LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -type ${j} -g -h -s ${data} > ${final_result_path}/${datasets[$i]}_hybrid_tree_${delta_type[j]}.txt
  done
done


##### code for different delta
# for ar in 0 0.25 0.5 0.75 1
# do
# data=/public/home/tangwei/graphdata/${dir_prefix[4]}/csr.${csr_suffix[4]}
# LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -r $ar -s ${data} > ${result_dir}/crossdelta/flickr-a${ar}.txt
# done

# for dr in 0.0001 0.0005 0.001 0.005 0.01 0.05
# do
# data=/public/home/tangwei/graphdata/${dir_prefix[4]}/csr.${csr_suffix[4]}
# LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ${exe} -n ${delta_num} -f $dr -s ${data} > ${result_dir}/crossdelta/flickr-d${dr}.txt
# done


