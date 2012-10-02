#!/usr/bin/env python
#coding: utf8

# Quick python program to plot chains produced by the
# save_chain() function in MC Julia.

from numpy import loadtxt
from matplotlib import pyplot as plot
from sys import argv, exit
from os.path import exists

if len(argv) < 2:
    print "Error: no filename given."
    print "Usage:   %s <filename>" % (argv[0])
    exit()
if len(argv) > 2:
    print "Warning: too many filenames."
    print "The following were ignored:  %s" % " ".join(argv[2:])

filename = argv[1]

if not exists(filename):
    print "Error: file not found: %s" % filename
    exit()

try:
    chain = loadtxt(filename, delimiter=",")
except ValueError:
    print "Error: could not read file %s." % filename
    exit()

s = chain.shape
if len(s) > 1:
    steps = s[1]
    dim = s[0]
else:
    steps = s[0]
    dim = 1

plot.figure(0)

# Loop over the chains for each parameter
for i in xrange(dim):
    # Plot chain
    plot.figure(0)
    sub = dim*100 + 10 + i+1
    plot.subplot(sub)
    if dim > 1:
        plot.plot(chain[i,:], color="black")
    else:
        plot.plot(chain[:], color="black")
    
    # Plot histogram
    plot.figure(i+1)
    plot.title("Histogram for parameter #%d" % (i+1))
    if dim > 1:
        plot.hist(chain[i,:], 30, color="black")
    else:
        plot.hist(chain[:], 30, color="black")

# Make correlation plot
if dim <= 9 and dim > 1:
    plot.figure(dim+1)
    for i in xrange(0, dim-1):
        for j in xrange(0, dim):
            if i < j:
                plot.subplot(dim-1, dim-1, i*(dim-1) + j)
                plot.ylabel("Parameter #%d" % (i+1))
                plot.xlabel("Parameter #%d" % (j+1))
                plot.scatter(chain[j,:], chain[i,:], 0.1, color="black")
plot.show()
