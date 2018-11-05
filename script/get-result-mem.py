#!/usr/bin/env python3

import os
import re

pathbase = "/public/home/tangwei/project/ligra/result/"
codes = ["full", "sharing", "chaining"]

def getfactor(filename):
  r = re.search(r"f(0\.\d+)", filename)
  return r.group(1)

def getmems(filename):
  ret = []
  with open(filename) as f:
    for line in f.readlines():
      if line.find("current") != -1 :
        ret.append(line.split()[-2])
  return ret

for code in [1]:
  rfile = pathbase + "vresult.csv"
  with open(rfile, "w") as rf:
    # directory = pathbase + "vt-" + code + "/"
    directory = pathbase 
    for filename in os.listdir(directory):
      if "vgraph" not in filename:
        continue
      factor = getfactor(filename)
      rf.write(factor + ",")
      rf.write(",".join(getmems(directory + filename)))
      rf.write("\n")
      
