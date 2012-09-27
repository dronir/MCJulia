#!/usr/bin/env python
#coding: utf8

from numpy import loadtxt
from matplotlib import pyplot as plot
from sys import argv, exit

filename = argv[1]

burn_in = 400

chain = loadtxt(filename, delimiter=",")
s = chain.shape
k = s[1]
dim = s[0]

for i in xrange(dim):
    plot.figure(0)
    sub = dim*100 + 10 + i+1
    plot.subplot(sub)
    plot.plot(chain[i,burn_in:], color="black")
    
    plot.figure(i+1)
    plot.hist(chain[i,burn_in:], 20, color="black")

plot.figure(0)
plot.axvline(burn_in, color="red")
plot.show()