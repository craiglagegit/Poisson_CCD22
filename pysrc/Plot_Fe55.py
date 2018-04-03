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
print outputfiledir+"/plots/Fe55_Sim_%d_%d.pdf"%(NumRuns,NumSpots)
dirbase = outputfiledir.split('run')[0]

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
Nx = ConfigData['PixelBoundaryNx']
Ny = ConfigData['PixelBoundaryNy']

Fe55RepulsionMult = ConfigData["Fe55RepulsionMult"]

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
        filename = dirbase+'run_%d/%s_%d_CC.dat'%(spot,outputfilebase,run)
        elec = ReadCCFile(filename, Nx, Ny)
        cfgfile = dirbase+'run_%d/pixel.cfg'%spot
        SpotConfigData = ReadConfigFile(cfgfile)
        xsum = 0.0
        ysum = 0.0
        elecsum = 0.0
        for i in range(Nx):
            for j in range(Ny):
                spotlist.data[i,j,spot] = elec[i,j]
                elecsum += elec[i,j]
                xsum = spotlist.x[i] * elec[i,j]
                ysum = spotlist.y[j] * elec[i,j]
        spotlist.xoffset[spot] = xsum / elecsum
        spotlist.yoffset[spot] = ysum / elecsum

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

# Now plot the results
# Create the output directory if it doesn't exist

xindex = Nx / 2
yindex = Ny / 2

AveCounts = spotlist.data.sum(axis=2) / float(NumSpots)
NormFactor = AveCounts[xindex, yindex] / ADU_correction
gauss = zeros([Nx, Ny])
print AveCounts.sum()

for i in range(Nx):
    for j in range(Ny):
        gauss[i, j] = Area(spotlist.x[i]-0.5, spotlist.x[i]+0.5, spotlist.y[j]-0.5, spotlist.y[j]+0.5, sigmaxs[0], sigmays[0], NormFactor)

if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

file = open(outputfiledir+"/fe55.txt","w")
file.write("Sigmax = %.4f, Sigmay = %.4f\n"%(sigmaxs[0], sigmays[0]))
file.close()
figure()
suptitle("Fe55 spot, %d spots, Fe55RepulsionMult = %.1f"%(NumSpots,Fe55RepulsionMult))
subplot(1,2,1)
title("X-Slice")
scatter(spotlist.x[:], AveCounts[:,yindex], marker='o', color='green', label = 'Sims')
scatter(spotlist.x[:], gauss[:,yindex], marker='x', color='red', label = 'Gaussian Fit')
xlabel('X (Pixels)')
ylabel('ADU')
ylim(0.0, 1.10 * AveCounts[xindex,yindex])
xlim(spotlist.xmin, spotlist.xmax)
legend(loc='center right')
text(-2.0,AveCounts[xindex,yindex] * 0.70, "SigmaX = %.4f"%sigmaxs[0], fontsize=16, color='blue')

subplot(1,2,2)
title("Y-Slice")
scatter(spotlist.y[:], AveCounts[xindex,:], marker='o', color='green', label='Sims')
scatter(spotlist.y[:], gauss[xindex,:], marker='x', color='red', label = 'Gaussian Fit')
xlabel('Y (Pixels)')
ylabel('ADU')
ylim(0.0, 1.10 * AveCounts[xindex,yindex])
xlim(spotlist.ymin, spotlist.ymax)
legend(loc='center right')
text(-2.0,AveCounts[xindex,yindex] * 0.70, "SigmaY = %.4f"%sigmays[0], fontsize=16, color='blue')
subplots_adjust(wspace=0.5)
savefig(outputfiledir+"/plots/Fe55_Sim_%d.pdf"%(NumSpots))


figure()
subplot(1,1,1,aspect = 1)
suptitle("Fe55 spot, %d spots, Fe55RepulsionMult = %.1f"%(NumSpots,Fe55RepulsionMult))

for i in range(Nx):
    plot([spotlist.x[i]-0.5, spotlist.x[i]-0.5], [spotlist.y[0]-0.5, spotlist.y[Ny-1]+0.5], color = 'black', ls = '--')
    plot([spotlist.x[i]+0.5, spotlist.x[i]+0.5], [spotlist.y[0]-0.5, spotlist.y[Ny-1]+0.5], color = 'black', ls = '--')
for j in range(Ny):
    plot([spotlist.x[0]-0.5, spotlist.x[Nx-1]+0.5], [spotlist.y[j]-0.5, spotlist.y[j]-0.5], color = 'black', ls = '--')
    plot([spotlist.x[0]-0.5, spotlist.x[Nx-1]+0.5], [spotlist.y[j]+0.5, spotlist.y[j]+0.5], color = 'black', ls = '--')


for i in range(Nx):
    for j in range(Ny):
        text(spotlist.x[i] - 0.30, spotlist.y[j] + 0.15, "%.1f"%AveCounts[i,j], color = 'green') 
        text(spotlist.x[i] - 0.30, spotlist.y[j] - 0.25, "%.1f"%gauss[i,j], color = 'red') 

text(spotlist.x[Nx-1] + 1.2, spotlist.y[yindex] + 0.15, "Sims", color = 'green') 
text(spotlist.x[Nx-1] + 1.2, spotlist.y[yindex] - 0.25, "Gaussian Fit", color = 'red') 
text(spotlist.x[Nx-1] + 1.2, spotlist.y[yindex] - 0.55, "$\sigma_x$=%.3f"%sigmaxs[0],color = 'red') 
text(spotlist.x[Nx-1] + 1.2, spotlist.y[yindex] - 0.85, "$\sigma_y$=%.3f"%sigmays[0],color = 'red') 

savefig(outputfiledir+"/plots/Fe55_AveCounts_%d.pdf"%NumSpots)

