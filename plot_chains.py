#!/usr/bin/env python
#coding: utf8

from numpy import loadtxt
from matplotlib import pyplot as plot
from sys import argv, exit

filename = argv[1]

burn_in = 0

chain = loadtxt(filename, delimiter=",")
s = chain.shape
if len(s) > 1:
    steps = s[1]
    dim = s[0]
else:
    steps = s[0]
    dim = 1


for i in xrange(dim):
    plot.figure(0)
    sub = dim*100 + 10 + i+1
    plot.subplot(sub)
    if dim > 1:
        plot.plot(chain[i,burn_in:], color="black")
    else:
        plot.plot(chain[burn_in:], color="black")
    
    plot.figure(i+1)
    if dim > 1:
        plot.hist(chain[i,burn_in:], 20, color="black")
    else:
        plot.hist(chain[burn_in:], 20, color="black")

plot.figure(0)
plot.axvline(burn_in, color="red")
plot.show()