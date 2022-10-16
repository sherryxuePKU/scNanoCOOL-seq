#-*- coding:utf-8 -*-
from __future__ import division
import sys
import numpy as np

depth = 10
infile = sys.argv[1]
f_infile = open(infile, "r")
for line in f_infile:
    line = line.strip()
    f = line.split()
    val_1 = np.array(f[5].split(","), dtype="int").max()
    val_2 = np.array(f[6].split(","), dtype="int").max()
    if val_1 + val_2 >= depth:
        print line

f_infile.close()
