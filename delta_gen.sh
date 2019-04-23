
# for i in 0 0.25 0.5 0.75 1 
# do
#   LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ./bin/BellmanFord -n 100 -r $i  -s /public/home/tangwei/graphdata/flickr/csr.flickr-growth >> tmp/flickr2.txt &
# done


for j in 0.0001
do
LD_PRELOAD=/public/home/tangwei/ytw/lib64/libstdc++.so.6 ./bin/BellmanFord -n 100 -s /public/home/tangwei/graphdata/twitter/csr.twitter_rv.net >> tmp/twitter.txt
done
