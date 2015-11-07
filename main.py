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
from functions import Showdata
from functions import Simulation
from functions import Compose
from functions import Decompose

#import genes
from genes import argP
from genes import cspD
from genes import diaA
from genes import dinB
from genes import dnaA
from genes import dnaB
from genes import dnaC
from genes import dnaE
from genes import dnaG
from genes import dnaN
from genes import dnaQ
from genes import dnaT
from genes import dnaX
from genes import gspB
from genes import gyrA
from genes import gyrB
from genes import hda
from genes import holA
from genes import holB
from genes import holC
from genes import holD
from genes import holE
from genes import polA
from genes import polB
from genes import priA
from genes import priB
from genes import priC
from genes import rarA
from genes import rdgC
from genes import recF
from genes import recG
from genes import recQ
from genes import rep
from genes import rnhA
from genes import rob
from genes import sbmC
from genes import seqA
from genes import topA
from genes import topB
from genes import tus
from genes import umuC
from genes import umuD
from genes import uvrD

#Sequence data
seq, mod, time, SubList = Reactions().Readseq('sequence.txt')

#Substance generate
Reactions().Monomer('YciV',60000,SubList)#RNA/ssDNA exonuclease 5'->3'specific(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6634-MONOMER)
Reactions().Monomer('YdaV',60000,SubList)#Rac prophage; predicted DNA replication protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6684-MONOMER)
Reactions().Monomer('YcdX',60000,SubList)#zinc-binding phosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6540-MONOMER)
Reactions().Monomer('CrfC',60000,SubList)#Clamp-binding sister replication fork colocalization protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11210-MONOMER)
Reactions().Monomer('RarA',60000,SubList)#recombination factor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12690-MONOMER)
Reactions().Monomer('Hda',60000,SubList)#regulator of DnaA that prevents premature reinitiation of DNA replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G7313-MONOMER)
Reactions().Monomer('DiaA',60000,SubList)#DnaA initiator-associating factor for replication initiation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=YRAO-MONOMER)
Reactions().Monomer('HolB',60000,SubList)#DNA polymerase III, delta prime subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11500-MONOMER)
Reactions().Monomer('HolA',60000,SubList)#DNA polymerase III, delta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11412-MONOMER)
Reactions().Monomer('Tus',60000,SubList)#DNA-binding protein; inhibition of replication at Ter sites(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11038-MONOMER)
Reactions().Monomer('DnaC',60000,SubList)#chromosome replication; initiation and chain elongation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10237-MONOMER)
Reactions().Monomer('Dam',60000,SubList)#DNA adenine methyltransferase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10204-MONOMER)
Reactions().Monomer('DnaG',60000,SubList)#DNA primase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10239-MONOMER)
Reactions().Monomer('DnaN',60000,SubList)#DNA polymerase III, beta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10242-MONOMER)
Reactions().Monomer('DnaT',60000,SubList)#primosomal protein DnaT(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10244-MONOMER)
Reactions().Monomer('Rep',60000,SubList)#Rep helicase, a single-stranded DNA dependent ATPase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10837-MONOMER)
Reactions().Monomer('PriB',60000,SubList)#primosomal replication protein N(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10764-MONOMER)
Reactions().Monomer('DnaQ',60000,SubList)#DNA polymerase III, epsilon subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10243-MONOMER)
Reactions().Monomer('DnaA',60000,SubList)#chromosomal replication initiator protein DnaA; DNA-binding transcriptional dual regulator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD03831)
Reactions().Monomer('MukE',60000,SubList)#protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11252-MONOMER)
Reactions().Monomer('MukB',60000,SubList)#cell division protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10618-MONOMER)
Reactions().Monomer('MukF',60000,SubList)#MukF dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12165-MONOMER)
Reactions().Monomer('NrdD',60000,SubList)#ribonucleoside-triphosphate reductase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=RIBONUCLEOSIDE-TRIP-REDUCT-MONOMER)
Reactions().Monomer('RecQ',60000,SubList)#ATP-dependent DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10833-MONOMER)
Reactions().Monomer('DinB',60000,SubList)#DNA polymerase IV (Y-family DNA polymerase; translesion DNA synthesis)(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6115-MONOMER)
Reactions().Monomer('NrdB',60000,SubList)#ribonucleoside diphosphate reductase 1, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDB-MONOMER)
Reactions().Monomer('NrdE',60000,SubList)#ribonucleoside-diphosphate reductase 2, alpha subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDE-MONOMER)
Reactions().Monomer('NrdF',60000,SubList)#ribonucleoside-diphosphate reductase 2, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDF-MONOMER)
Reactions().Monomer('HolC',60000,SubList)#DNA polymerase III, chi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11413-MONOMER)
Reactions().Monomer('HolD',60000,SubList)#DNA polymerase III, psi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11414-MONOMER)
Reactions().Monomer('LigB',60000,SubList)#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11334-MONOMER)
Reactions().Monomer('PriA',60000,SubList)#primosome factor N'(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10763-MONOMER)
Reactions().Monomer('LigA',60000,SubList)#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10534-MONOMER)
Reactions().Monomer('SbcC',60000,SubList)#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
Reactions().Monomer('DnaE',60000,SubList)#DNA polymerase III, alpha subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
Reactions().Monomer('Ssb',60000,SubList)#ssDNA-binding protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10976-MONOMER)
Reactions().Monomer('SbcD',60000,SubList)#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11094-MONOMER)
Reactions().Monomer('DnaB',60000,SubList)#replicative DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10236-MONOMER)
Reactions().Monomer('HolE',60000,SubList)#DNA polymerase III, theta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11505-MONOMER)
Reactions().Monomer('RecF',60000,SubList)#ssDNA and dsDNA binding, ATP binding(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10828-MONOMER)
Reactions().Monomer('LexA',60000,SubList)#LexA DNA-binding transcriptional repressor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD00205)
Reactions().Monomer('DnaK',60000,SubList)#chaperone protein DnaK(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10241-MONOMER)
Reactions().Monomer('DnaJ',60000,SubList)#chaperone protein DnaJ(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10240-MONOMER)
Reactions().Monomer('MutT',60000,SubList)#8-oxo-dGTP diphosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10626-MONOMER)
Reactions().Monomer('DnaX',60000,SubList)#DNA polymerase III, tau subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10245-MONOMER)
Reactions().Monomer('GspB',60000,SubList)#calcium-binding protein required for initiation of chromosome replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11263-MONOMER)
Reactions().Monomer('UvrD',60000,SubList)#ssDNA translocase and dsDNA helicase - DNA helicase II(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11064-MONOMER)
Reactions().Monomer('PolA',60000,SubList)#DNA polymerase I, 5'-->3'polymerase, 5'-->3'and3'-->5'exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10746-MONOMER)
Reactions().Monomer('NrdA',60000,SubList)#ribonucleoside diphosphate reductase 1, alpha atom fun subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDA-MONOMER)
Reactions().Complex('DNA polymerase III, core enzyme',0,SubList)
Reactions().Complex('DNA polymerase III, preinitiation complex',0,SubList)
Reactions().Complex('DNA polymerase III, beta subunit',0,SubList)
Reactions().Complex('DNA polymerase III, tau subunit dimer',0,SubList)
Reactions().Complex('DNA polymerase III, psi-chi subunit',0,SubList)
Reactions().Enzyme('DNA polymerase III, holoenzyme',0,SubList)

#Gillespie Test
logt, logd, t, tend = Showdata().logger(time, SubList, 0, 0.001)
events = [Compose('DNA polymerase III, core enzyme',['DnaE','DnaQ','HolE'],[1,1,1],0.001),
          #Compose('DNA polymerase III, preinitiation complex',['DnaX','HolB','HolA'],[3,1,1],0.1),
          Compose('DNA polymerase III, beta subunit',['DnaN'],[2],0.001),
          #Compose('DNA polymerase III, tau subunit dimer',['DnaX'],[2],0.1),
          #Compose('DNA polymerase III, psi-chi subunit',['HolC','HolD'],[1,1],0.1),
          Decompose('DNA polymerase III, core enzyme',['DnaE','DnaQ','HolE'],[1,1,1],1),
          #Decompose('DNA polymerase III, preinitiation complex',['DnaX','HolB','HolA'],[3,1,1],5),
          Decompose('DNA polymerase III, beta subunit',['DnaN'],[2],1)]
          #Decompose('DNA polymerase III, tau subunit dimer',['DnaX'],[2],5),
          #Decompose('DNA polymerase III, psi-chi subunit',['HolC','HolD'],[1,1],5)]
Simulation().run(t, tend, SubList, events, logt, logd)
Showdata().figure("DnaE", logt, logd, SubList)
Showdata().figure("DnaQ", logt, logd, SubList)
Showdata().figure("HolE", logt, logd, SubList)
#Showdata().figure("DnaX", logt, logd, SubList)
#Showdata().figure("HolB", logt, logd, SubList)
#Showdata().figure("HolA", logt, logd, SubList)
Showdata().figure("DnaN", logt, logd, SubList)
#Showdata().figure("HolC", logt, logd, SubList)
#Showdata().figure("HolD", logt, logd, SubList)
Showdata().figure("DNA polymerase III, core enzyme", logt, logd, SubList)
Showdata().figure("DNA polymerase III, beta subunit", logt, logd, SubList)
Showdata().save("complex.png")
"""
#Polymerization Process Test
r = [0]
i = 0
time = []
er = []
process = len(seq)
for location in range(process):
    dnaB().execute(location,mod,SubList)
for location in range(process):
    holB().execute(location,mod,SubList,r)
    time.append(i)
    er.append(r[0])
    i += 1

#Reading Chain5'-3'

#Ragging Chain5'-3'

#Save result
Simulation().Wcplot(time,er)
Simulation().Wcwrite(mod)
"""
Simulation().Makedata()
