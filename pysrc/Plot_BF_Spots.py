#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 10-Sep17

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import os, sys, time, subprocess, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines
from scipy.special import erf
from scipy.optimize import fmin_powell
from scipy import stats
sys.path.append(os.path.realpath('/global/homes/c/cslage/Software/forward_model_varying_i'))
import forward

#****************SUBROUTINES*****************
def FOM(params):
    global spotlist
    [sigmax, sigmay] = params
    result = forward.forward(spotlist,sigmax,sigmay)
    return result

#****************MAIN PROGRAM*****************
global spotlist
# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
NumSpots = int(sys.argv[3])
NumRuns = int(sys.argv[4])

ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
print outputfiledir
print outputfiledir+"/plots/BF_Sim_%d_%d.pdf"%(NumRuns,NumSpots)
dirbase = outputfiledir.split('bfrun')[0]

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
Nx = ConfigData['PixelBoundaryNx']
Ny = ConfigData['PixelBoundaryNy']
stampxmin = -(int(Nx/2)+0.5)
stampxmax = -stampxmin
stampymin = -(int(Ny/2)+0.5)
stampymax = -stampymin

imaxs = zeros([NumRuns])
sigmaxs = zeros([NumRuns])
sigmays = zeros([NumRuns])

for run in range(NumRuns):
    spotlist = Array2dSet(stampxmin,stampxmax,Nx,stampymin,stampymax,Ny,NumSpots)
    for spot in range(NumSpots):
        filename = dirbase+'bfrun_%d/%s_%d_CC.dat'%(spot,outputfilebase,run)
        elec = ReadCCFile(filename, Nx, Ny)
        cfgfile = dirbase+'bfrun_%d/bf.cfg'%spot
        SpotConfigData = ReadConfigFile(cfgfile)
        spotlist.xoffset[spot] = SpotConfigData['Xoffset'] / SpotConfigData['PixelSizeX']
        spotlist.yoffset[spot] = SpotConfigData['Yoffset'] / SpotConfigData['PixelSizeY']
        for i in range(Nx):
            for j in range(Ny):
                spotlist.data[i,j,spot] = elec[i,j]

    # Now run the forward modeling
    param0 = [1.00, 1.00]
    args = ()
    Result = fmin_powell(FOM, param0, args)
    print Result
    imax = spotlist.imax.mean()
    ADU_correction = Area(-0.5,0.5,-0.5,0.5,Result[0],Result[1],1.0)
    sigmaxs[run] = Result[0]
    sigmays[run] = Result[1]
    imaxs[run] = imax * ADU_correction
    del(spotlist)


# Now plot the results
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

file = open(outputfiledir+"/bf.txt","w")
figure()
title("Baseline - Sigmax = Sigmay = 1.0, Offsets=random, With Diffusion")
scatter(imaxs, sigmaxs, color = 'green', lw = 2, label = 'Sigma-x')
scatter(imaxs, sigmays, color = 'red', lw = 2, label = 'Sigma-y')

slope, intercept, r_value, p_value, std_err = stats.linregress(imaxs[10:-1],sigmaxs[10:-1])
xplot=linspace(0.0,150000.0,100)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='blue', lw = 2, ls = '--')
txslope = slope * 100.0 * 50000.0
file.write("X Slope = %.2f %% per 50K e-, Intercept = %.3f\n"%(txslope,intercept))
slope, intercept, r_value, p_value, std_err = stats.linregress(imaxs[10:-1],sigmays[10:-1])
xplot=linspace(0.0,150000.0,100)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='black', lw = 2, ls = '--')
tyslope = slope * 100.0 * 50000.0
file.write("Y Slope = %.2f %% per 50K e-, Intercept = %.3f\n"%(tyslope,intercept))
xlabel('Central Peak(electrons)')
ylabel('Sigma (Pixels)')
legend(loc= 'lower right')
ylim(0.95, 1.10)
xlim(0.0,150000.0)
xticks([0.0,50000,100000])
xslope_meas = 0.71
xslope_sigma = 0.09
yslope_meas = 0.94
yslope_sigma = 0.07
chi_squared = (((txslope - xslope_meas) / xslope_sigma)**2 + ((tyslope - yslope_meas) / yslope_sigma)**2) / 2.0
text(10000.0,1.08,"%d Simulated spots"%NumSpots, fontsize=16, color='blue')
text(10000.0,1.07,"X Slope = %.2f %% per 50K e-, Data = %.2f +/- %.2f"%(txslope, xslope_meas, xslope_sigma), fontsize=16, color='blue')
text(10000.0,1.06,"Y Slope = %.2f %% per 50K e-, Data = %.2f +/- %.2f"%(tyslope, yslope_meas, yslope_sigma), fontsize=16, color='blue')
text(10000.0,1.05,"$\chi^2$ = %.2f"%chi_squared, fontsize=16, color='blue')

savefig(outputfiledir+"/plots/BF_Sim_%d_%d.pdf"%(NumRuns,NumSpots))
file.close()

