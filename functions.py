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
import glob

#General reactions
class Reactions:
    def __init__(self,time=0):
        self.time = time
        self.sl = []
        self.evl = []

    def setup(self):
        print "------------------------------------------------------------------"
        print "| "+RED+"A whole-genome based simulation of prokaryotic DNA replication"+ENDC+" |"
        print "| "+GREEN+"MODEL : "+YELLOW+"Escherichia coli K-12 substr. MG1655                  "+ENDC+" |"
        print "------------------------------------------------------------------"
        return self.time, self.sl, self.evl

    def Readseq(self, readseq, sublist):
        seq = open(readseq)
        self.f = seq.read()
        self.m = np.zeros((len(self.f),len(sublist)+5))
        for i in range(len(self.f)):
            if self.f[i] == 'a':
                self.m[i][len(sublist)+2] = 0.23
                self.m[i][len(sublist)+4] = 0.73
            if self.f[i] == 't':
                self.m[i][len(sublist)+2] = 0.24
                self.m[i][len(sublist)+4] = 0.74
            if self.f[i] == 'g':
                self.m[i][len(sublist)+2] = 0.25
                self.m[i][len(sublist)+4] = 0.75
            if self.f[i] == 'c':
                self.m[i][len(sublist)+2] = 0.26
                self.m[i][len(sublist)+4] = 0.76
            if self.f[i] == 'M':
                self.m[i][len(sublist)+2] = 0.88
                self.m[i][len(sublist)+4] = 0.88
        return self.f,self.m

    def Monomer(self,name,num,sublist):
        self.name = str(name)
        self.num = num
        print "{:<20}".format(YELLOW+self.name+ENDC)+" = "+BLUE+`self.num`+ENDC
        gn = [self.name,self.num]
        sublist.append(gn)
        return sublist

    def Complex(self,name,num,complexlist):
        self.name = str(name)
        self.num = num
        print "{:<55}".format(RED+self.name+ENDC)+" = "+BLUE+`self.num`+ENDC
        gn = [self.name,self.num]
        complexlist.append(gn)
        return complexlist

    def Region(self,region,state):
        pass

    def Getindex(self, name, sublist):
        for i in range(len(sublist)):
            if sublist[i][0] == name:
                return i
        if name == "ds":
            return len(sublist)

    def Increase(self, location, name, mod, sublist):
        mod[location][Reactions().Getindex(name,sublist)] += 1
        return mod

    def Decrease(self, location, name, mod, sublist):
        mod[location][Reactions().Getindex(name,sublist)] -= 1
        return mod

    def Events(self, ev, evlist):
        evlist.append(ev)
        return evlist

class Compose:
    def __init__(self,name,components,comnum,k):
        self.name = name
        self.components = components
        self.comnum = comnum
        self.k = k

    def propensity(self, state):
        self.p = self.k
        for i in range(len(self.components)):
            self.p = self.p * state[Reactions().Getindex(self.components[i],state)][1] ** self.comnum[i]
        return self.p

    def execute(self, state):
        for i in range(len(self.components)):
            state[Reactions().Getindex(self.components[i],state)][1] -= self.comnum[i]
        state[Reactions().Getindex(self.name,state)][1] += 1
        return state

class Decompose:
    def __init__(self,name,components,comnum,k):
        self.name = name
        self.components = components
        self.comnum = comnum
        self.k = k

    def propensity(self, state):
        self.p = self.k * state[Reactions().Getindex(self.name,state)][1]
        return self.p

    def execute(self, state):
        for i in range(len(self.components)):
            state[Reactions().Getindex(self.components[i],state)][1] += self.comnum[i]
        state[Reactions().Getindex(self.name,state)][1] -= 1
        return state

class Bind:
    def __init__(self,sub,sublist,location,k):
        pass

    def propensity(self, state, location, k):
        return

    def execute(self, state, location, k):
        return

class Unbind:
    def __init__(self,sub,sublist,location,k):
        pass

    def propensity(self, state, location, k):
        return

    def execute(self, state, location, k):
        return

class Showdata:
    def __init__(self):
        pass

    def logger(self, t, sublist, st, et):
        otime = [t]
        data = []
        for i in range(len(sublist)):
            data.append([sublist[i][1]])
        self.data = data
        self.time = otime
        self.st = st
        self.et = et
        return self.time, self.data, self.st, self.et

    def getdata(self, t, sublist, logt, logd):
        for i in range(len(sublist)):
            logd[i].append(sublist[i][1])
        logt.append(t)
        return logt, logd

    def figure(self, name, logt, logd, sublist):
        for i in range(len(sublist)):
            if sublist[i][0] == str(name):
                plt.plot(logt, logd[i][:], label=sublist[i][0])

    def save(self, name="figure.png"):
        plt.legend()
        plt.savefig(name)
        plt.close()

    def png(self,figlist,logt,logd,sublist,pngname):
        for name in figlist:
            Showdata().figure(name, logt, logd, sublist)
        if pngname == 'default':self.pngname = figlist[-1]+'.png'
        else: self.pngname = pngname+'.png'
        Showdata().save(self.pngname)

class Simulation:
    def __init__(self):
        pass

    def step(self, time, state, events):
        atotal = 0
        alist = []
        for i in range(len(events)):
                atotal += events[i].propensity(state)
                alist.append(events[i].propensity(state))
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

    def run(self, t, tend, SubList, events, logt, logd, mod):
        while t <= tend:
            t, SubList = Simulation().step(t, SubList, events)
            logt, logd = Showdata().getdata(t, SubList, logt, logd)
        return t, SubList, logt, logd

    def Wcplot(self, time, data):
        plt.plot(time, data, 'b', label="ERROR BASE")
        plt.xlabel("Time[s]", fontsize=12)
        plt.ylabel("ERROR BASE", fontsize=12)
        plt.title("Error Accumulation", fontsize=14)
        plt.savefig("error.png")
        plt.close()

    def Save(self, mod):
        np.save('result.npy', mod)
        f = open('loadresult.py','w')
        f.write("#!/usr/bin/env python\n")
        f.write("import numpy as np\n")
        f.write("f = np.load('result.npy')\n")
        f.write("print f\n")
        f.close()

    def Makedata(self, dirname):
        if dirname == 'default':dirname = 'result'
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
        if os.path.exists(pwd+'/result.npy'):shutil.move('result.npy',pwd+"/"+dirname)
        if os.path.exists(pwd+'/loadresult.py'):shutil.move('loadresult.py',pwd+"/"+dirname)
        for name in glob.glob('*.png'):
            shutil.move(name,pwd+"/"+dirname)
