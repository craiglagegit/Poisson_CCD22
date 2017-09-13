#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 26-Jan-15

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
from matplotlib.colors import from_levels_and_colors
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
    
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Vpl = ConfigData["Vparallel_lo"]
Vph = ConfigData["Vparallel_hi"]
# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]

ZMult = 2.0
kmax = int(dat.Elec.shape[2]*0.90)
nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1
nxcenter = nxx/2
nycenter = nyy/2

dnx = ScaleFactor*GridsPerPixelX/2
dny = ScaleFactor*GridsPerPixelY/2
csxcenter = nxcenter + dnx
bgycenter = nycenter + dny

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

ChargeDepth(dat, outputfiledir+'/charge.txt', 88*ScaleFactor, 45*ScaleFactor, dnx, dny/2, kmax)
carriers = ['Electron', 'Hole', 'Mobile', 'Fixed']
plotdatas = [dat.Elec, dat.Hole, dat.Hole-dat.Elec, dat.rho]
cmap0 = cm.get_cmap("jet")
cmap1 = cm.get_cmap("seismic")
cmaps = [cmap0, cmap0, cmap1, cmap1]
ForceZeros = [False, False, True, True]
xslicemins = [nxcenter-dnx, nxcenter-dnx, nxcenter-dnx, nxcenter-dnx]
xslicemaxs = [nxcenter+dnx, nxcenter+dnx, nxcenter+dnx, nxcenter+dnx]

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

for i, plotdata in enumerate(plotdatas):

    fig = figure(figsize = (12,12))
    suptitle("%s Charge Distribution"%carriers[i], fontsize = 36)

    ax1=axes([0.10,0.40,0.50,0.50],aspect=1)
    ax1.set_title("X-Y Slice")
    ax1.set_xticks([])
    ax1.set_yticks([])
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 2, nxmin, nxmax, nymin, nymax, 0, kmax, ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax1.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax2=axes([0.10,0.20,0.50,0.20], aspect=1)
    ax2.set_title("X-Z Slice")
    ax2.set_xticks([])
    ax2.set_yticks([0.0,2.0,4.0])
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 1, nxmin, nxmax, 0, kmax, nymin, nymax, ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax2.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax3=axes([0.10,0.10,0.50,0.10])
    ax3.set_title("X-Cut", fontsize=24)
    ax3.set_xlabel("X (Microns)", fontsize=18)
    ax3.set_ylabel("Charge Density \n(arb. units)", fontsize=18)
    ax3.set_yticks([0])
    [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 0, nxmin, nxmax, nycenter-1, nycenter+2, 0, kmax)
    ax3.plot(xpoints, plotslice)
    
    ax4=axes([0.60,0.40,0.20,0.50], aspect=1)
    ax4.set_title("Y-Z Slice")
    ax4.set_xticks([0.0,2.0,4.0])
    ax4.set_yticks([])
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 0, 0, kmax, nymin, nymax, xslicemins[i], xslicemaxs[i], ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax4.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')        
    ax4.invert_xaxis()

    ax5=axes([0.80,0.40,0.10,0.50])
    ax5.set_title("Y-Cut", fontsize=24)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    for tick in ax5.get_xticklabels():
        tick.set_rotation(90)
    ax5.set_ylabel("Y (Microns)", fontsize=18)
    ax5.set_xlabel("Charge Density \n(arb. units)", fontsize=18)
    ax5.set_xticks([0])
    [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 1, nymin, nymax, nxcenter-1, nxcenter+2, 0, kmax)
    ax5.plot(plotslice, xpoints)
    ax5.set_xlim(ax5.get_xlim()[::-1])

    ax6=axes([0.65,0.20,0.10,0.10])
    ax6.set_title("Z-Cut", fontsize=24)
    ax6.set_ylabel("Log Charge \n Density \n(arb. units)", fontsize=18)
    ax6.set_xticks([0.0,1.0,2.0])
    ax6.set_xlabel("Z ( Microns)", fontsize=18)
    ax6.yaxis.set_label_position("right")
    ax6.set_ylim(-4.0, 1.0)
    ax6.set_ylim(-4.0, 1.0)
    ax6.text(3.0,-7.0,"Z-Axis has a %.1fX Scale Multiplier"%ZMult, fontsize = 18)
    if i == 2:
        [plotslicen, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicep, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, csxcenter-1, csxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-8.0,"+ Charge in Red, - Charge in Blue", fontsize = 18)        
    elif i == 3:
        [plotslicep, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicen, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, csxcenter-1, csxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-8.0,"+ Charge in Red, - Charge in Blue", fontsize = 18)        
    else:
        [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxmin, nxmax, nymin, nymax)
        ax6.plot(xpoints, log10(plotslice/plotslice.max()+1.0E-5))
    ax6.set_xlim(ax6.get_xlim()[::-1])
    ax6.text(3.0,-9.0,"Oxide in yellow, Gates in green", fontsize = 18)        
    savefig(outputfiledir+"/plots/%sDistribution_New_%d.pdf"%(carriers[i],run))
    close(fig)
