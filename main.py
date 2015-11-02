#!/usr/bin/env python
#color setting
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

#import public module
import sys
import math
import re
import csv
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *

#import whole-ecoli replication module
from functions import Reactions
from functions import Enzyme
from functions import Propensity
from functions import Simulation

print RED+"This is whole-ecoli replication"+ENDC

#Sequence data
seq,mod,time = Reactions().Readseq('sequence.txt')
print GREEN+"time"+ENDC+" = "+BLUE+`time`+ENDC

#Substance setting
SubList = []                                                #ID
SubList.append(Reactions().Generate('dnaA',100))            #0
SubList.append(Reactions().Generate('dnaB',100))            #1
SubList.append(Reactions().Generate('dnaC',200))            #2
SubList.append(Reactions().Generate('dnaG',0))              #3
SubList.append(Reactions().Generate('RNaseH',0))            #4
SubList.append(Reactions().Generate('SSB',0))               #5
SubList.append(Reactions().Generate('Topo1',0))             #6
SubList.append(Reactions().Generate('SDC',0))               #7
SubList.append(Reactions().Generate('DNApol1',0))           #8
SubList.append(Reactions().Generate('DNApol3',100))         #9
SubList.append(Reactions().Generate('DNApol3holoenzyme',0)) #10
SubList.append(Reactions().Generate('OriC9',5))             #11
SubList.append(Reactions().Generate('OriC13',3))            #12
SubList.append(Reactions().Complex('dnaB','dnaC',0))        #13
SubList.append(Reactions().Complex('dnaA','dnaB/dnaC',0))   #14
SubList.append(Reactions().Complex('dnaA','ATP',0))         #15
SubList.append(Reactions().Complex('OriC9','dnaA',0))       #16
SubList.append(Reactions().Complex('OriC13','dnaA',0))      #17

#for i in range(50):
#    SubList = Reactions().Compose(SubList,'dnaA','dnaB',0.05)
#for i in range(50):
#    SubList = Reactions().Decompose(SubList,'dnaA/dnaB',5.00)

#for i in mod:
#    print i
#print
r = [0]
i = 0
time = []
er = []
process = len(seq)
for location in range(process):
    Enzyme().dnaB(location,mod,SubList,1)
for location in range(process):
    Enzyme().DNApol3(location,mod,SubList,1,r)
    time.append(i)
    er.append(r[0])
    i += 1

#Save result
plt.plot(time,er,'b',label="ERROR BASE")
plt.xlabel("Time[s]", fontsize=12)
plt.ylabel("ERROR BASE", fontsize=12)
plt.title("Error Accumulation", fontsize=14)
plt.savefig("error.png")
result = open('result.txt', 'w')
for line in mod:
    result.write(str(line)+'\n')
result.close()
Simulation().Makedata()
