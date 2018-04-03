#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 10-Sep17

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
skip = int(sys.argv[3])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
NumElec = ConfigData["NumElec"]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_Pts.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print "No data in Pts file.  Quitting"
    sys.exit()
lines.remove(lines[0])

if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the electron paths
    print "Making array electron path plots\n"
    # Plotting the paths along a line through the center

    NumSteps = 2000
    xpaths = zeros([NumElec, NumSteps])
    zxpaths = zeros([NumElec, NumSteps])
    ypaths = zeros([NumElec, NumSteps])
    zypaths = zeros([NumElec, NumSteps])
    maxsteps = zeros([NumElec], dtype=int)
    
    vertical_zoom = 0.5
    figure()
    suptitle("Electron Path Plot - Vertical Zoom = %.1f"%vertical_zoom, fontsize = 24)
    subplots_adjust(wspace=0.2)

    for line in lines:
        values = line.split()
        id = int(values[0])
        step = int(values[1])
        phase = int(values[2])
        if step > NumSteps - 1:
            continue
        
        x = float(values[3])
        y = float(values[4])
        z = float(values[5])
        if isnan(x) or isnan(y) or isnan(z):
            continue
        xpaths[id, step] = x
        zxpaths[id, step] = z            
        ypaths[id, step] = y
        zypaths[id, step] = z
        maxsteps[id] = step
    for id in range(0, NumElec, skip):
        subplot(1,2,1,aspect=vertical_zoom)
        plot(xpaths[id,0:maxsteps[id]], zxpaths[id,0:maxsteps[id]], color='k', linewidth = 0.01)
        subplot(1,2,2,aspect=vertical_zoom)
        plot(ypaths[id,0:maxsteps[id]], zypaths[id,0:maxsteps[id]], color='k', linewidth = 0.01)

    subplot(1,2,1,aspect=vertical_zoom)
    ylabel("Z(microns)")
    xlabel("X (microns)")
    xticks([20,30,40,50])
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0]+PixelSizeX, ConfigData["PixelBoundaryUpperRight"][0]-PixelSizeX)
    subplot(1,2,2,aspect=vertical_zoom)
    ylabel("Z(microns)")
    xlabel("Y (microns)")
    xticks([20,30,40,50])
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][1]+PixelSizeY, ConfigData["PixelBoundaryUpperRight"][1]-PixelSizeY)
    savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths_%d.pdf"%run)

