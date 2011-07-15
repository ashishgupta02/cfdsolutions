#!/usr/bin/env python

# ../License.h
from numpy import *

# If the package has been installed correctly, this should work:
import Gnuplot, Gnuplot.funcutils
import os, sys
from optparse import OptionParser

#-----------------------------------------------------
def VariableID(var):
     if (var.lower() == "rho"):
	return 2
     if (var.lower() == "rhou"):
	return 3
     if (var.lower() == "rhov"):
	return 4
     if (var.lower() == "rhow"):
	return 5
     if (var.lower() == "rhoet"):
	return 6
     if (var.lower() == "res"):
	return 7

#-----------------------------------------------------
def ResidualPlot(g, filename, variable):
    length = len(variable)
    if length <= 3:
        iDim = 1
        jDim = length
    else:
        iDim = 2
        jDim = 3
    # Get the Modified Time Stamp of the file
    g('set style function lines')
    g('set size 1.0, 1.0')
    g('set origin 0.0, 0.0')
    g('set multiplot')
    for p in range(length):
        size1 = 1.0/float(iDim)
        size2 = 1.0/float(jDim)        
        string = 'set size ' + str(size1) + ', ' + str(size2)
        g(string)
        if p < 3:
            origin1 = 0.0
            origin2 = p*(1.0/float(jDim))
        else:
            origin1 = 1.0/float(iDim)
            origin2 = (p-3)*(1.0/float(jDim))  
        string = 'set origin ' + str(origin1) + ', ' + str(origin2)
        g(string)
        g('set grid')
        g('unset key')
        string = variable[p] + ' Residual Plot'
        g.title(string)
        g.xlabel('Iterations')
        string = variable[p]
        g.ylabel(string)
        g('set logscale y')
        string = '"' + filename + '" using 1:' + str(VariableID(variable[p])) + ' ls ' + str(p+1)
        g.plot(string)
    g('unset multiplot')

#-----------------------------------------------------
def LoopPlot(filename, variable):
    # Get the Modified Time Stamp of the file
    mtime1 = os.stat(filename).st_mtime;    
    g = Gnuplot.Gnuplot(debug=0)
    ResidualPlot(g, filename, variable)
    while True:
        mtime2 = os.stat(filename).st_mtime
    	if mtime1 != mtime2:
            ResidualPlot(g, filename, variable)
	    mtime1 = mtime2
    raw_input('Please press return to continue...\n')
    g.reset()

#-----------------------------------------------------
def main():
    usage = "usage: %prog -f Filename -v Variable -p"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", action="store", type="string", dest="filename", default="residual.res")
    parser.add_option("-v", "--var",  action="append", type="string", dest="variable")
    (options, args) = parser.parse_args()

    if options.filename == None:
        print "ERROR: Input File Not Provided !"
        sys.exit()

    if options.variable == None:
        print "ERROR: No variable Supplied to Parse !"
        sys.exit()

    # Check if File exits
    if not os.path.exists(options.filename):
        print "ERROR: Input File Does Not Exists !"
        sys.exit()
    
    # Check if Parameter file exists
    LoopPlot(options.filename, options.variable)
    return

#-----------------------------------------------------
# Main:
if __name__ == '__main__':
    main()

