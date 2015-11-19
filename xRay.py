#!/usr/bin/env python

BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'
from numpy.random import *

f = open('sequence.txt','r')
target = f.read()
mutationseq = ''
STRENGTH = input(GREEN+"X-ray STRENGTH (%)"+ENDC+" = ")#%
MUTATION = 0
for base in target:
    er = rand() * 100
    if STRENGTH >= er:
        mutationseq += 'M'
        MUTATION += 1
    else:
        mutationseq += base
f.close()
f2 = open('mutationseq.txt','w')
f2.write(mutationseq)
f2.close()
print RED + 'MUTATION '+ENDC+ '= ' + `float(MUTATION/float(len(target)))*100` + GREEN + ' %' + ENDC
