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
        self.db = []
        self.original = open('sequence.txt','r').read()

    def setup(self):
        print "------------------------------------------------------------------"
        print "| "+RED+"A whole-genome based simulation of prokaryotic DNA replication"+ENDC+" |"
        print "| "+GREEN+"MODEL : "+YELLOW+"Escherichia coli K-12 substr. MG1655                  "+ENDC+" |"
        print "------------------------------------------------------------------"
        return self.time, self.sl, self.evl

    def Target(self,tg):
        self.tg = tg
        return self.tg

    def Initdata(self):
        return self.db

    def addDB(self,gene,protein,resion,list):
        list.append([gene,protein,Reactions().Getbase(resion)])
        return list

    def Readseq(self, readseq, sublist):
        seq = open(readseq)
        self.f = seq.read()
        self.m = np.zeros((len(self.f),len(sublist)+5))
        for i in range(len(self.f)):
            if self.f[i] == 'a':
                self.m[i][len(sublist)+1] = 0.23#a(r)
                self.m[i][len(sublist)+3] = 0.76#t(l)
            if self.f[i] == 't':
                self.m[i][len(sublist)+1] = 0.26#t(r)
                self.m[i][len(sublist)+3] = 0.73#a(l)
            if self.f[i] == 'g':
                self.m[i][len(sublist)+1] = 0.24#g(r)
                self.m[i][len(sublist)+3] = 0.75#c(l)
            if self.f[i] == 'c':
                self.m[i][len(sublist)+1] = 0.25#c(r)
                self.m[i][len(sublist)+3] = 0.74#g(l)
            if self.f[i] == 'M':
                self.m[i][len(sublist)+1] = 0.28#mutation(r)
                self.m[i][len(sublist)+3] = 0.78#mutation(l)
        return self.f,self.m

    def Baseinfo(self, ID, mod):
            self.m = mod
            if round(ID,0) == 0.0: self.rl = 'r'
            if round(ID,0) == 1.0: self.rl = 'l'
            number,base = [3,4,5,6,8],['a','g','c','t','M']
            for i in range(len(number)):
                if round((ID-0.2)*100,0) == number[i] or round((ID-0.7)*100,0) == number[i]: self.base = base[i]
            return self.rl,self.base
            #print Reactions().Baseinfo(mod[i][Reactions().Getindex('r',SubList)],mod)[1]
            #print Reactions().Baseinfo(0.23,mod)

    def Monomer(self,name,num,sublist,target):
        self.name = str(name)
        self.jd,self.gene = Database(target).Match(self.name)
        self.num = num
        if self.jd == 1:
            print GREEN+'OK '+ENDC+': '+GREEN+"{:<5}".format(self.gene)+ENDC+"{:<3}".format('->')+GREEN+"{:<5}".format(self.name)+ENDC+"{:>33}".format(' = ')+BLUE+`self.num`+ENDC
        else:
            self.num = 0
            print RED+'NG '+ENDC+': '+RED+"{:<5}".format(self.gene)+ENDC+"{:<3}".format('-|')+RED+"{:<5}".format(self.name)+ENDC+"{:>33}".format(' = ')+BLUE+`self.num`+ENDC
        gn = [self.name,self.num]
        sublist.append(gn)
        return sublist

    def Complex(self,name,num,complexlist):
        self.name = str(name)
        self.num = num
        print GREEN+'OK '+ENDC+': '+"{:<52}".format(YELLOW+self.name+ENDC)+' = '+BLUE+`self.num`+ENDC
        gn = [self.name,self.num]
        complexlist.append(gn)
        return complexlist

    def Region(self,region,state):
        pass

    def Getindex(self, name, sublist):
        for i in range(len(sublist)):
            if sublist[i][0] == name: return i
        if name == "ds": return len(sublist)
        if name == "r": return len(sublist)+1
        if name == "l": return len(sublist)+3

    def Getbase(self, rg):
        return self.original[rg[0]-1:rg[1]-1]

    def Increase(self, location, name, mod, sublist):
        mod[location][Reactions().Getindex(name,sublist)] += 1
        return mod

    def Decrease(self, location, name, mod, sublist):
        mod[location][Reactions().Getindex(name,sublist)] -= 1
        return mod

    def Events(self, ev, evlist):
        evlist.append(ev)
        return evlist

class Database:
    def __init__(self,target):
        self.gcellDB = Reactions().Initdata()
        Reactions().addDB('yciV','YciV',[1321244,1322125],self.gcellDB)
        Reactions().addDB('ydaV','YdaV',[1420007,1420753],self.gcellDB)
        Reactions().addDB('ycdX','YcdX',[1098102,1098839],self.gcellDB)
        Reactions().addDB('yjcZ','CrfC',[4327383,4328261],self.gcellDB)
        Reactions().addDB('rarA','RarA',[937217,938560],self.gcellDB)
        Reactions().addDB('hda','Hda',[2616097,2616798],self.gcellDB)
        Reactions().addDB('diaA','DiaA',[3293831,3294421],self.gcellDB)
        Reactions().addDB('holB','HolB',[1154985,1155989],self.gcellDB)
        Reactions().addDB('holA','HolA',[669797,670828],self.gcellDB)
        Reactions().addDB('tus','Tus',[1682283,1683212],self.gcellDB)
        Reactions().addDB('dnaC','DnaC',[4598261,4598998],self.gcellDB)
        Reactions().addDB('dam','Dam',[3513099,3513935],self.gcellDB)
        Reactions().addDB('dnaG','DnaG',[3209129,3210874],self.gcellDB)
        Reactions().addDB('dnaN','DnaN',[3879244,3880344],self.gcellDB)
        Reactions().addDB('dnaT','DnaT',[4599001,4599540],self.gcellDB)
        Reactions().addDB('rep','Rep',[3958700,3960721],self.gcellDB)
        Reactions().addDB('priB','PriB',[4423543,4423857],self.gcellDB)
        Reactions().addDB('dnaQ','DnaQ',[236067,236798],self.gcellDB)
        Reactions().addDB('dnaA','DnaA',[3880349,3881752],self.gcellDB)
        Reactions().addDB('mukE','MukE',[974845,975549],self.gcellDB)
        Reactions().addDB('mukB','MukB',[975549,980009],self.gcellDB)
        Reactions().addDB('mukF','MukF',[973542,974864],self.gcellDB)
        Reactions().addDB('nrdD','NrdD',[4458545,4460683],self.gcellDB)
        Reactions().addDB('recQ','RecQ',[4003887,4005716],self.gcellDB)
        Reactions().addDB('dinB','DinB',[250898,251953],self.gcellDB)
        Reactions().addDB('nrdB','NrdB',[2345406,2346536],self.gcellDB)
        Reactions().addDB('nrdE','NrdE',[2799370,2801514],self.gcellDB)
        Reactions().addDB('nrdF','NrdF',[2801524,2802483],self.gcellDB)
        Reactions().addDB('holC','HolC',[4481860,4482303],self.gcellDB)
        Reactions().addDB('holD','HolD',[4605826,4606239],self.gcellDB)
        Reactions().addDB('ligB','LigB',[3817511,3819193],self.gcellDB)
        Reactions().addDB('priA','PriA',[4122635,4124833],self.gcellDB)
        Reactions().addDB('ligA','LigA',[2526183,2528198],self.gcellDB)
        Reactions().addDB('sbcC','SbcC',[411831,414977],self.gcellDB)
        Reactions().addDB('dnaE','DnaE',[205126,208608],self.gcellDB)
        Reactions().addDB('ssB','Ssb',[4272148,4272684],self.gcellDB)
        Reactions().addDB('sbcD','SbcD',[414974,416176],self.gcellDB)
        Reactions().addDB('dnaB','DnaB',[4262337,4263752],self.gcellDB)
        Reactions().addDB('holE','HolE',[1923132,1923362],self.gcellDB)
        Reactions().addDB('recF','RecF',[3878171,3879244],self.gcellDB)
        Reactions().addDB('lexA','LexA',[4255138,4255746],self.gcellDB)
        Reactions().addDB('dnaK','DnaK',[12163,14079],self.gcellDB)
        Reactions().addDB('dnaJ','DnaJ',[14168,15298],self.gcellDB)
        Reactions().addDB('mutT','MutT',[111044,111433],self.gcellDB)
        Reactions().addDB('dnaX','DnaX',[491316,493247],self.gcellDB)
        Reactions().addDB('gspB','GspB',[3451530,3451949],self.gcellDB)
        Reactions().addDB('uvrD','UvrD',[3996006,3998168],self.gcellDB)
        Reactions().addDB('polA','PolA',[4044989,4047775],self.gcellDB)
        Reactions().addDB('nrdA','NrdA',[2342887,2345172],self.gcellDB)
        Reactions().addDB('priC','PriC',[489509,490036],self.gcellDB)
        self.seq = open(target,'r').read()

    def Match(self, protein):
        for i in self.gcellDB:
            if protein == i[1]:
                if i[2] in self.seq: return 1,i[0]
                else: return 0,i[0]

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
        if atotal == 0:
            print RED+'\ncan\'t continue '+ENDC+': '+BLUE+'Inviable environment'+ENDC
            sys.exit()
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

    def Run(self, t, tend, SubList, events, logt, logd, mod):
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

    def Save(self, mod, sublist):
        np.save('result.npy', mod)
        label = []
        for name in sublist:
            label.append(name[0])
        label.append('Double/Single Strand')
        label.append('ReadingOriginal')
        label.append('ReadingReplicated')
        label.append('LaggingOriginal')
        label.append('LaggingReplicated')
        label = ",".join(label)
        f = open('npy2csv.py','w')
        f.write("#!/usr/bin/env python\n")
        f.write("import numpy as np\n")
        f.write("f = np.load('result.npy')\n")
        f.write("np.savetxt('result.csv', f, fmt='%.02f',delimiter=',',header=\""+label+"\")\n")
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
                if os.path.exists(pwd+"/"+dirname):dirname = raw_input(RED+"ERROR "+GREEN+"Please input other name : "+ENDC)
                else: break
        os.mkdir(dirname)
        os.system('chmod +x npy2csv.py')
        if os.path.exists(pwd+'/result.npy'):shutil.move('result.npy',pwd+"/"+dirname)
        if os.path.exists(pwd+'/npy2csv.py'):shutil.move('npy2csv.py',pwd+"/"+dirname)
        for name in glob.glob('*.png'):shutil.move(name,pwd+"/"+dirname)
