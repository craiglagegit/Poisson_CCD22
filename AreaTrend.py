#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 16-Nov-15

#This program plots the pixel area plots from the Poisson CCD solver
import matplotlib
matplotlib.use("Agg")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************

def ReadConfigFile(filename):
    # This reads the config file for the necessary settings
    # and returns a dictionary with the values
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    ConfigData = {}
    try:
        for line in lines:
            ThisLine=line.strip().split()
            ThisLineLength=len(ThisLine)
            if ThisLineLength < 3:
                continue
            if list(ThisLine[0])[0]=='#' or ThisLine[0]=='\n':
                continue
            try:
                ParamName = ThisLine[0]
                ThisLine.remove(ThisLine[0])
                for counter,item in enumerate(ThisLine):
                    if list(item)[0] == '#':
                        del ThisLine[counter:] # Strip the rest of the line as a comment
                        continue
                    if item == '=':
                        ThisLine.remove(item)
                        continue
                if len(ThisLine) == 0:
                    continue
                elif len(ThisLine) == 1:
                    ThisParam = ThisLine[0]
                    try: ConfigData[ParamName] = int(ThisParam)
                    except ValueError:
                        try: ConfigData[ParamName] = float(ThisParam)
                        except ValueError:
                            try:
                                ConfigData[ParamName] = ThisParam
                            except ValueError:
                                print "Error reading .cfg file"
                else:
                    ThisParam = []
                    for item in ThisLine:
                        try: ThisParam.append(int(item))
                        except ValueError:
                            try: ThisParam.append(float(item))
                            except ValueError:
                                ThisParam.append(item)
                    ConfigData[ParamName] = ThisParam
            except:
                continue
    except:
        print "Error reading .cfg file"

    return ConfigData

def ReadAreaFile(filename, nx, ny):
    area = zeros([nx, ny])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        items = line.split()
        if items[0] == "Nx":
            continue
        i = int(items[0])
        j = int(items[1])
        area[i,j] = float(items[2])
    return area


#****************MAIN PROGRAM*****************
dirs = ["data/run1","data/run2","data/run3","data/run4","data/run6"]
Area_0 = 100.0
central_areas = []
plus_x_areas = []
plus_y_areas = []
plus_x_plus_y_areas = []
charges = []

for dir in dirs:
    configfile = dir+"/bf.cfg"
    ConfigData = ReadConfigFile(configfile)
    outputfilebase = ConfigData["outputfilebase"]
    outputfiledir = ConfigData["outputfiledir"]
    Nx = ConfigData["PixelBoundaryNx"]
    Ny = ConfigData["PixelBoundaryNy"]
    NxCenter = int((ConfigData["FilledPixelCoords_0_0"][0] - ConfigData["PixelBoundaryLowerLeft"][0]) / ConfigData["PixelSize"])
    NyCenter = int((ConfigData["FilledPixelCoords_0_0"][1] - ConfigData["PixelBoundaryLowerLeft"][1]) / ConfigData["PixelSize"])
    NumElec = ConfigData["CollectedCharge_0_0"]

    filename = outputfiledir + '/' + outputfilebase +'_0_Area' + '.dat'

    area = ReadAreaFile(filename, Nx, Ny)

    central_areas.append(area[4,4] - Area_0)
    plus_x_areas.append(area[5,4] - Area_0)
    plus_y_areas.append(area[4,5] - Area_0)
    plus_x_plus_y_areas.append(area[5,5] - Area_0)
    charges.append(float(NumElec)/1000)

figure()
subplots_adjust(hspace=0.5, wspace = 0.5)
suptitle("Pixel Areas vs Charge")
subplot(2,2,3)
title("Pixel (0,0)")
scatter(charges, central_areas)
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(charges,central_areas)
xplot=linspace(0.0, 200.0, 2.0)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(10.0, -25.0, "Slope = %.4f"%slope)
text(10.0, -0.8*25.0, "R_value = %.4f"%r_value)
xlabel("Charge(Ke-)")
ylabel("$\delta Area (\mu ^2)$")
xlim(0.0,200.0)
xticks([0,100,200])
ylim(0.0, -30.0)

subplot(2,2,4)
title("Pixel (1,0)")
scatter(charges, plus_x_areas)
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(charges,plus_x_areas)
xplot=linspace(0.0, 200.0, 2.0)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(10.0, -0.4, "Slope = %.4f"%slope)
text(10.0, -0.8*0.4, "R_value = %.4f"%r_value)
xlabel("Charge(Ke-)")
ylabel("$\delta Area (\mu ^2)$")
xlim(0.0,200.0)
xticks([0,100,200])
ylim(0.0, -0.5)

subplot(2,2,1)
title("Pixel (0,1)")
scatter(charges, plus_y_areas)
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(charges,plus_y_areas)
xplot=linspace(0.0, 200.0, 2.0)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(10.0, 2.5, "Slope = %.4f"%slope)
text(10.0, 0.8*2.5, "R_value = %.4f"%r_value)
xlabel("Charge(Ke-)")
ylabel("$\delta Area (\mu ^2)$")
xlim(0.0,200.0)
xticks([0,100,200])
ylim(0.0, 3.0)

subplot(2,2,2)
title("Pixel (1,1)")
scatter(charges, plus_x_plus_y_areas)
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(charges,plus_x_plus_y_areas)
xplot=linspace(0.0, 200.0, 2.0)
yplot = slope * xplot + intercept
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(10.0, 0.8, "Slope = %.4f"%slope)
text(10.0, 0.8*0.8, "R_value = %.4f"%r_value)
xlabel("Charge(Ke-)")
ylabel("$\delta Area (\mu ^2)$")
xlim(0.0,200.0)
xticks([0,100,200])
ylim(0.0, 1.0)

savefig("Area_Trend.png")
