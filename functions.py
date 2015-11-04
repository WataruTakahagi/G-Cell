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

#General reactions
class Reactions:
    def __init__(self,time=0):
        self.time = time

    def Readseq(self, readseq='sequence.txt'):
        seq = open(readseq)
        self.f = seq.read()
        self.m = []
        for i in range(len(self.f)):
            self.m.append([0,0,0,0,0,0,0,0,0,0,self.f[i],''])
        #ds,dnaA,dnaB,dnaC,dnaG,SSB,Topo1,SDC,DNApol1,DNApol3
        return self.f,self.m,self.time

    def Generate(self,name,num):
        self.name = str(name)
        self.num = num
        print "{:<27}".format(YELLOW+self.name+ENDC)+" = "+BLUE+`self.num`+ENDC
        gn = [self.name,self.num]
        return gn

    def Complex(self,name1,name2,num):
        self.name1 = str(name1)
        self.name2 = str(name2)
        self.Cname = str(self.name1+"/"+self.name2)
        self.num = num
        print "{:<27}".format(RED+self.Cname+ENDC)+" = "+BLUE+`self.num`+ENDC
        return [self.Cname,self.num]

class Compose:
    def __init__(self, csubA, csubB):
        self.csubA = csubA
        self.csubB = csubB

    def propensity(self, state, k):
        self.k = k
        for i in range(len(state)):
            if state[i][0] == str(self.csubA): self.a = state[i][1]
            if state[i][0] == str(self.csubB): self.b = state[i][1]
        self.p = self.k * self.a * self.b
        return self.p

    def execute(self, state):
        self.Complex = str(self.csubA+'/'+self.csubB)
        for i in range(len(state)):
            if state[i][0] == str(self.csubA):
                state[i][1] = state[i][1] - 1
                print '['+GREEN+str(state[i][0])+ENDC+','+BLUE+str(state[i][1])+ENDC+']',
            if state[i][0] == str(self.csubB):
                state[i][1] = state[i][1] - 1
                print '['+GREEN+str(state[i][0])+ENDC+','+BLUE+str(state[i][1])+ENDC+']',
            if state[i][0] == self.Complex:
                state[i][1] = state[i][1] + 1
                print '['+GREEN+str(state[i][0])+ENDC+','+RED+str(state[i][1])+ENDC+']',
        print ENDC
        return state

class Decompose:
    def __init__(self, dsubAB):
        self.dsubAB = dsubAB
        self.Complex = dsubAB.split('/')

    def propensity(self, state, k):
        self.k = k
        for i in range(len(state)):
            if state[i][0] == str(self.dsubAB): self.ab = state[i][1]
        self.p = self.k * self.ab
        return self.p

    def execute(self, state):
        for i in range(len(state)):
            if state[i][0] == str(self.Complex[0]):
                state[i][1] = state[i][1] + 1
                print "["+GREEN+str(state[i][0])+ENDC+","+RED+str(state[i][1])+ENDC+"]",
            if state[i][0] == str(self.Complex[1]):
                state[i][1] = state[i][1] + 1
                print "["+GREEN+str(state[i][0])+ENDC+","+RED+str(state[i][1])+ENDC+"]",
            if state[i][0] == self.dsubAB:
                state[i][1] = state[i][1] - 1
                print "["+GREEN+str(state[i][0])+ENDC+","+BLUE+str(state[i][1])+ENDC+"]",
        print ENDC
        return state

class dnaA:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class dnaB:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        mseq[location][0] = 1
        state[1][1] -= 1
        return location, mseq, state

class dnaC:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class dnaG:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class RNaseH:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class SSB:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class Topo1:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class SDC:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        mseq[location][7] = 1
        state[7][1] -= 1
        return location, mseq, state

class DNApol1:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class DNApol3:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state, r):
        if mseq[location][0] == 1:
            mseq[location][9] = 1
            state[9][1] -= 1
            error = 0.001#ERROR RATE(%)
            rate = rand()*100
            if mseq[location][10]=='a':
                if rate >= error:mseq[location][11]='t'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='t':
                if rate >= error:mseq[location][11]='a'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='g':
                if rate >= error:mseq[location][11]='c'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='c':
                if rate >= error:mseq[location][11]='g'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            mseq[location][0] = 0
            location += 1
        return location, mseq, state, r

class DNApol3holoenzyme:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, state):
        return

class Simulation:
    def __init__(self):
        pass

    def Step(self, time, state, events, k):
        atotal = 0
        alist = []
        for i in range(len(events)):
                atotal += events[i].propensity(state, k[i])
                alist.append(events[i].propensity(state, k[i]))
        tau = float((1/atotal)*math.log1p(1/rand()))
        newt = time + tau
        a0, l = 0, 0
        r = rand()*atotal
        for a in alist:
            a0 += a
            if a0 > r: j = l
            else: l += 1
        news = events[j].execute(state)
        return newt, news

    def Wcplot(self, time, data):
        plt.plot(time, data, 'b', label="ERROR BASE")
        plt.xlabel("Time[s]", fontsize=12)
        plt.ylabel("ERROR BASE", fontsize=12)
        plt.title("Error Accumulation", fontsize=14)
        plt.savefig("error.png")

    def Wcwrite(self, mod):
        result = open('result.txt', 'w')
        for line in mod:
            result.write(str(line)+'\n')
        result.close()

    def Makedata(self, dirname="result"):
        pwd = os.getcwd()
        if os.path.exists(pwd+"/"+dirname):
            swt = 1
            print BLUE+dirname+RED+" already exists !!"+ENDC
            dirname = raw_input(YELLOW+"Please input other name : "+ENDC)
            while swt == 1:
                if os.path.exists(pwd+"/"+dirname):
                    dirname = raw_input(RED+"ERROR "+GREEN+"Please input other name : "+ENDC)
                else: break
        os.mkdir(dirname)
        if os.path.exists(pwd+'/error.png'): shutil.move('error.png',pwd+"/"+dirname)
        if os.path.exists(pwd+'/result.txt'):shutil.move('result.txt',pwd+"/"+dirname)
