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

data = [0.77, 0.91, 0.005943, 0.01227, -1.275, [0.001591, 0.002966, 0.000994, 0.000882, 0.000862, 0.000501]]
dataij = [0, 0, 1, 1, 0, [4,2,5,4,5,8]]

# First, read the .cfg file
run = 0

configfile = sys.argv[1]
step = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
NxCenter = 4
NyCenter = 4
Area_0 = 100.0

filename = outputfiledir + '/' + outputfilebase +'_%d_Area'%step+'.dat'

area = ReadAreaFile(filename, Nx, Ny)

NumElec = step * 1000

x = []
y = []
xfit = []
yfit = []
xneg = []
yneg = []
fullfile = open(outputfiledir+"/full_corr.txt","w")
for i in range(Nx):
    for j in range(Ny):
        rsquared = float((i - NxCenter)**2 + (j - NyCenter)**2)
        yvalue = (area[i,j] - Area_0) / Area_0
        if (i == NxCenter and j == NyCenter) or rsquared > 18.0:
            continue
        if yvalue < 0:
            xneg.append(rsquared)
            yneg.append(-yvalue)
        else:
            x.append(rsquared)
            y.append(yvalue)
        if rsquared > 1.1 and abs(i - NxCenter) < 3 and abs(j - NyCenter) < 3 and yvalue > 0:
            fullfile.write("i = %d, j = %d, C = %.6f\n"%(i-NxCenter,j-NxCenter,yvalue))
            xfit.append(rsquared)
            yfit.append(yvalue)

xdata = []
ydata = []
for i,dat in enumerate(data):
    if i != 5:
        if dataij[i] == 0:
            continue
        else:
            xdata.append(dataij[i])
            ydata.append(data[i])
    else:
        for j,d in enumerate(dat):
            xdata.append(dataij[i][j])
            ydata.append(data[i][j])

figure()
title("Covariance vs Distance: %d e-"%NumElec)
xscale('log')
yscale('log')
xlim(0.8, 100.0)
ylim(1E-5, 1E-1)
scatter(xdata, ydata, marker = 'x', s = 50, color = 'green', label = 'Data')
scatter(array(x), array(y), color = 'blue', label = 'Sim')
#scatter(xfit, yfit, marker="x", color = 'black')
#scatter(xneg, yneg, color = 'magenta', label = 'Negative Values')

from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(log10(xfit),log10(yfit))
xplot=linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
text(10.0, 0.005, "Slope = %.3f"%slope)
text(10.0, 0.00316, "C10 = %.5g"%((area[NxCenter+1,NyCenter] - Area_0) / Area_0))
text(10.0, 0.002, "C01 = %.5g"%((area[NxCenter,NyCenter+1] - Area_0) / Area_0))
xlabel("$i^2 + j^2$")
ylabel("$\delta$ Area / Area")
xlim(0.8,100)
legend()
savefig(outputfiledir+"/plots/Area_UnBiased_Corr_%d_%d.pdf"%(run, step))
#savefig("Area_%d_%d.pdf"%(run, step))

file = open(outputfiledir+"/corr.txt","w")
file.write("Slope = %.3f\n"%slope)
file.write("C10 = %.5g\n"%((area[NxCenter+1,NyCenter] - Area_0) / Area_0))
file.write("C01 = %.5g\n"%((area[NxCenter,NyCenter+1] - Area_0) / Area_0))
fullfile.write("Slope = %.3f\n"%slope)
fullfile.write("C10 = %.5g\n"%((area[NxCenter+1,NyCenter] - Area_0) / Area_0))
fullfile.write("C01 = %.5g\n"%((area[NxCenter,NyCenter+1] - Area_0) / Area_0))
file.close()
fullfile.close()
figure()
title("Pixel Area: %d e-"%NumElec)

for i in range(Nx+1):
    plot([10.0 + 10.0 * i, 10.0 + 10.0 * i], [10.0, 100.0], color = 'black')
for j in range(Ny+1):
    plot([10.0, 100.0], [10.0 + 10.0 * j, 10.0 + 10.0 * j], color = 'black')
for i in range(Nx):
    for j in range(Ny):
        if i == NxCenter and j == NyCenter:
            textcolor = 'red'
        else:
            textcolor = 'black'
        text(11.0 + 10.0 * i, 14.0 + 10.0 * j, str(area[i,j]), color = textcolor, fontsize = 8, fontweight = 'bold')
xlim(10.0, 100.0)
ylim(10.0, 100.0)
xlabel("X(microns)")
ylabel("Y(microns)")
savefig(outputfiledir+"/plots/PixelAreas_UnBiased_Corr_%d_%d.pdf"%(run, step))
#savefig("PixelAreas_%d_%d.pdf"%(run, step))
