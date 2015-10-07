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
import numpy as np
import matplotlib.pyplot as plt

#import whole-ecoli replication module
from functions import Reactions
from functions import Enzyme
from functions import Simulation

print RED+"This is whole-ecoli replication"+ENDC

#Sequence data
seq,mod,time = Reactions().Readseq('sequence.txt')
print GREEN+"time"+ENDC+" = "+BLUE+`time`+ENDC
print mod
#Substance setting
SubList = []
SubList.append(Reactions().Generate('dnaA',100))
SubList.append(Reactions().Generate('dnaB',100))
SubList.append(Reactions().Generate('dnaC',200))
SubList.append(Reactions().Generate('dnaG',0))
SubList.append(Reactions().Generate('RNaseH',0))
SubList.append(Reactions().Generate('SSB',0))
SubList.append(Reactions().Generate('Topo1',0))
SubList.append(Reactions().Generate('SDC',0))
SubList.append(Reactions().Generate('DNApol1',0))
SubList.append(Reactions().Generate('DNApol3',0))
SubList.append(Reactions().Generate('DNApol3holoenzyme',0))
SubList.append(Reactions().Complex('dnaA','dnaB',20))
SubList.append(Reactions().Complex('dnaA','ATP',20))

for i in range(50):
    SubList = Reactions().Compose(SubList,'dnaA','dnaB',0.05)
for i in range(50):
    SubList = Reactions().Decompose(SubList,'dnaA/dnaB',5.00)

#Initiation
#Reactions().Compose(SubList,'dnaA','dnaB',0.05)

#Enzyme setting
#Enzyme().dnaA(0,1,100)
#Enzyme().dnaB(0,1,200)
#Enzyme().dnaC(0,1,300)
#Enzyme().dnaG(0,1,400)
#Enzyme().RNaseH(0,1,100)
#Enzyme().SSB(0,1,100)
#Enzyme().Topo1(0,1,100)
#Enzyme().SDC(0,1,100)
#Enzyme().DNApol1(0,1,200)
#Enzyme().DNApol3(0,1,100)
#Enzyme().DNApol3holoenzyme(0,1,100)

#Initiation
Enzyme().dnaA(0,1,100)
