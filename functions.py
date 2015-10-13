#!/usr/bin/env python

import sys
import math
import re
from numpy.random import *
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

class Reactions():
    def __init__(self,time=0):
        self.time = time

    def Readseq(self, readseq):
        seq = open(readseq)
        self.f = seq.read()
        self.m = []
        for i in range(len(self.f)):
            self.m.append([0,0,0,0,0])
        #self.m = re.sub(r'[a-z]',list(00),self.f)
        #self.m = list(self.m)
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

    def Compose(self, list, csubA, csubB, k):
        self.Complex = str(csubA+"/"+csubB)
        self.k = k
        for i in range(len(list)):
            if list[i][0] == str(csubA):
                self.a = list[i][1]
                list[i][1] = list[i][1] - 1
                print "["+GREEN+str(list[i][0])+ENDC+","+BLUE+str(list[i][1])+ENDC+"]",
            if list[i][0] == str(csubB):
                self.b = list[i][1]
                list[i][1] = list[i][1] - 1
                print "["+GREEN+str(list[i][0])+ENDC+","+BLUE+str(list[i][1])+ENDC+"]",
            if list[i][0] == self.Complex:
                list[i][1] = list[i][1] + 1
                print "["+GREEN+str(list[i][0])+ENDC+","+RED+str(list[i][1])+ENDC+"]",
        print ENDC
        self.p = self.k * self.a * self.b
        return list

    def Decompose(self, list, dsubAB,k):
        self.Complex = dsubAB.split('/')
        self.k = k
        for i in range(len(list)):
            if list[i][0] == str(self.Complex[0]):
                list[i][1] = list[i][1] + 1
                print "["+GREEN+str(list[i][0])+ENDC+","+RED+str(list[i][1])+ENDC+"]",
            if list[i][0] == str(self.Complex[1]):
                list[i][1] = list[i][1] + 1
                print "["+GREEN+str(list[i][0])+ENDC+","+RED+str(list[i][1])+ENDC+"]",
            if list[i][0] == dsubAB:
                self.ab = list[i][1]
                list[i][1] = list[i][1] - 1
                print "["+GREEN+str(list[i][0])+ENDC+","+BLUE+str(list[i][1])+ENDC+"]",
        print ENDC
        self.p = self.k * self.ab
        return list

class Enzyme():
    def __init__(self):
        pass

    def dnaA(self, time, state, k):
        return

    def dnaB(self, location, mseq, state, k):
        mseq[location][0] = 1
        state[1][1] -= 1
        return location, mseq, state

    def dnaC(self, time, state, k):
        return

    def dnaG(self, time, state, k):
        return

    def RNaseH(self, time, state, k):
        return

    def SSB(self, time, state, k):
        return

    def Topo1(self, time, state, k):
        return

    def SDC(self, time, state, k):
        return

    def DNApol1(self, time, state, k):
        return

    def DNApol3(self, time, state, k):
        return

    def DNApol3holoenzyme(self, time, state, k):
        return

class Simulation:
    def __init__(self):
        pass

    def Propensity(self, state, k):
        pass

    def Step(self, time, state, events):
        atotal = 0
        alist = []
        for i in range(len(events)):
            atotal += events[i]
