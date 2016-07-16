#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 16-Nov-15

#This program plots the pixel area plots from the Poisson CCD solver
import matplotlib
matplotlib.use("Agg")
from pylab import *
import sys, time, h5py
from scipy import stats, polyfit, poly1d

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

def ReadVertexFile(filename, nx, ny, NumAngles):
    vx = zeros([nx, ny, NumAngles])
    vy = zeros([nx, ny, NumAngles])
    theta = zeros([nx, ny, NumAngles])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    lasti = 3457
    lastj = 3457
    for line in lines:
        items = line.split()
        if items[0] == "X0":
            continue
        i = (int(float(items[0])) - int(1.5*ConfigData["PixelSize"])) / int(ConfigData["PixelSize"])
        j = (int(float(items[1])) - int(1.5*ConfigData["PixelSize"])) / int(ConfigData["PixelSize"])
        if i != lasti or j != lastj:
            k = 0
            lasti = i
            lastj = j

        #print i,j,k
        theta[i,j,k] = float(items[2])
        vx[i,j,k] = float(items[3])
        vy[i,j,k] = float(items[4])
        k += 1
    return (theta, vx, vy)

#****************MAIN PROGRAM*****************
dirs = ("data/run%d"%i for i in range(1,9))
Vx0 = 60.0
Vy0 = 60.0

Center_UR_xs = []
Center_UR_ys = []
Center_R_xs = []
Center_R_ys = []

zs = []

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
    NumAngles = 4 * ConfigData["NumVertices"] + 4

    filename = outputfiledir + '/' + outputfilebase +'_0_Vertices' + '.dat'

    (thetas, vxs, vys) = ReadVertexFile(filename, Nx, Ny, NumAngles)

    StartedNegatives = False
    for i, theta in enumerate(thetas[4,4,:]):
        if theta > 0:
            if StartedNegatives:
                Center_UR_xs.append(Vx0 - vxs[4,4,i+32])
                Center_UR_ys.append(Vy0 - vys[4,4,i+32])
                Center_R_xs.append(Vx0 - vxs[4,4,i+4])
                Center_R_ys.append(Vy0 - 5.0 - vys[4,4,i+4])
                #Center_R_ys.append(Vy0 - 5.0 + 5.0 * tan(9.0 * pi / 256.0) - vys[4,4,i+4])
                break
            else:
                continue
        else:
            StartedNegatives = True
    """
    StartedNegatives = False
    for i, theta in enumerate(thetas[4,5,:]):
        if theta > 0:
            if StartedNegatives:
                #if dir == "data/run6":
                #    print theta
                #    print vxs[5,4,i+32]
                #    print vys[5,4,i+32]
                #    sys.exit()
                Plus_Y_UR_xs.append(Vx0 - vxs[4,5,i+32])
                Plus_Y_UR_ys.append(Vy0 + 10.0 - vys[4,5,i+32])
                break
            else:
                continue
        else:
            StartedNegatives = True
    """
    zs.append(ConfigData["ElectronZ0"])

zs = array(zs)

figure()
subplots_adjust(hspace=0.5, wspace = 0.5)
suptitle("Vertex Shift vs Z0")
subplot(2,2,1)
title("Pixel (0,0) UR X Shift")
scatter(zs, Center_UR_xs)
xplot=linspace(3.0, 100.0, 100)
xfit = 12.0
yoff = 0.0
yplot = tanh(xplot/xfit)*Center_UR_xs[0] + yoff
predict = tanh(zs/xfit)*Center_UR_xs[0] + yoff
error = 0.0
for i, dat in enumerate(Center_UR_xs):
    error += abs((dat - predict[i]) / dat)
error = (error / float(len(Center_UR_xs))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.2, "tanh(z / %.2f)"%xfit)
text(20.0, 0.1, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
ylim(0.0, 0.5)

subplot(2,2,2)
title("Pixel (0,0) UR Y Shift")
scatter(zs, Center_UR_ys)
xplot=linspace(3.0, 100.0, 100)
xfit = 7.0
yoff = 0.2
yplot = tanh(xplot/xfit)*Center_UR_ys[0] + yoff
predict = tanh(zs/xfit)*Center_UR_ys[0] + yoff
error = 0.0
for i, dat in enumerate(Center_UR_ys):
    error += abs((dat - predict[i]) / dat)
error = (error / float(len(Center_UR_ys))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.3, "tanh(z / %.2f)"%xfit)
text(20.0, 0.1, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
#ylim(0.0, -0.4)

subplot(2,2,3)
title("Pixel (0,0) R X Shift")
scatter(zs, Center_R_xs/Center_R_xs[0])
xplot=linspace(3.0, 100.0, 100)
xfit = 10.0
yplot = tanh(xplot/xfit)
predict = tanh(zs/xfit)
error = 0.0
for i, dat in enumerate(Center_R_xs):
    error += abs((dat - predict[i] * Center_R_xs[0]) / dat)
error = (error / float(len(Center_R_xs))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.3, "tanh(z / %.2f)"%xfit)
text(20.0, 0.1, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
ylim(0.0, 1.2)

subplot(2,2,4)
title("Pixel (0,0) R Y Shift")
scatter(zs, Center_R_ys/Center_R_ys[0])
xplot=linspace(3.0, 100.0, 100)
xfit = 15.0
yplot = tanh(xplot/xfit)
predict = tanh(zs/xfit)
error = 0.0
for i, dat in enumerate(Center_R_ys):
    error += abs((dat - predict[i] * Center_R_ys[0]) / dat)
error = (error / float(len(Center_R_ys))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.3, "tanh(z / %.2f)"%xfit)
text(20.0, 0.1, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
ylim(0.0, 1.2)
"""
subplot(2,2,2)
title("Pixel (0,0) UR Y Shift")
scatter(zs, Center_UR_ys)
coeffs=polyfit(1.0/zs, Center_UR_ys/Center_UR_ys[0], 2)
print "(0,0) UR Y"
print coeffs
#coeffs = array([ 46.0, -18.9,   -1.28,    1.03])
#coeffs = array([-580.0, 410.0, -90.0, 2.69, 0.979])
polynom=poly1d(coeffs)
xplot=linspace(3.0, 100.0, 100)
yplot = polynom(1.0/xplot) * Center_UR_ys[0]
predict = polynom(1.0/zs) * Center_UR_ys[0]
error = 0.0
for i, dat in enumerate(Center_UR_ys):
    error += abs((dat - predict[i]) / dat)
error = (error / float(len(Center_UR_ys))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.5, "4th order poly fit in 1/z")
text(20.0, 0.45, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
#ylim(0.0, 1.0)

subplot(2,2,3)
title("Pixel (0,0) R X Shift")
scatter(zs, Center_R_xs)
coeffs=polyfit(1.0/zs, Center_R_xs/Center_R_xs[0], 2)
print "(0,0) R X"
print coeffs
#coeffs = array([ 46.0, -18.9,   -1.28,    1.03])
#coeffs = array([-580.0, 410.0, -90.0, 2.69, 0.979])
polynom=poly1d(coeffs)
xplot=linspace(3.0, 100.0, 100)
yplot = polynom(1.0/xplot) * Center_R_xs[0]
predict = polynom(1.0/zs) * Center_R_xs[0]
error = 0.0
for i, dat in enumerate(Center_R_xs):
    error += abs((dat - predict[i]) / dat)
error = (error / float(len(Center_R_xs))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.2, "4th order poly fit in 1/z")
text(20.0, 0.12, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
#ylim(0.0, 0.8)

subplot(2,2,4)
title("Pixel (0,0) R Y Shift")
scatter(zs, Center_R_ys)
coeffs=polyfit(1.0/zs, Center_R_ys/Center_R_ys[0], 2)
print "(0,0) R Y"
print coeffs
#coeffs = array([ 46.0, -18.9,   -1.28,    1.03])
#coeffs = array([-580.0, 410.0, -90.0, 2.69, 0.979])
polynom=poly1d(coeffs)
xplot=linspace(3.0, 100.0, 100)
yplot = polynom(1.0/xplot) * Center_R_ys[0]
predict = polynom(1.0/zs) * Center_R_ys[0]
error = 0.0
for i, dat in enumerate(Center_R_ys):
    error += abs((dat - predict[i]) / dat)
error = (error / float(len(Center_R_ys))) * 100
plot(xplot, yplot, color='red', lw = 2, ls = '--')
text(20.0, 0.0, "4th order poly fit in 1/z")
text(20.0, -0.01, "Ave error = %.2f %%"%error)
xlabel("Z(microns)")
ylabel("$Vertex Shift (\mu)$")
xlim(0.0,100.0)
xticks([0,50,100])
#ylim(0.0, 0.15)
"""

savefig("Vertex_Trend.png")
