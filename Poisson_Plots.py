#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Jan-15

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************
class Array3dHDF5(object):
    def __init__(self, dir, filebase, LogEField, run):
        phifile = dir+'/'+filebase+'_'+str(run)+'_phi'
        rhofile = dir+'/'+filebase+'_'+str(run)+'_rho'
        hdfphi = h5py.File(phifile,'r')
        Dimension = hdfphi[hdfphi.items()[0][0]].attrs[u'Dimension']
        self.nx=Dimension[0]
        self.ny=Dimension[1]
        self.nz=Dimension[2]
        
        Lower_Left = hdfphi[hdfphi.items()[0][0]].attrs[u'Lower_Left']
        self.xmin=Lower_Left[0]
        self.ymin=Lower_Left[1]
        self.zmin=Lower_Left[2]

        Upper_Right = hdfphi[hdfphi.items()[0][0]].attrs[u'Upper_Right']
        self.xmax=Upper_Right[0]
        self.ymax=Upper_Right[1]
        self.zmax=Upper_Right[2]
        
        self.dx=(self.xmax-self.xmin)/self.nx
        self.dy=(self.ymax-self.ymin)/self.ny
        self.dz=(self.zmax-self.zmin)/self.nz
        self.volume = self.dx * self.dy * self.dz
        
        self.x=linspace(self.xmin+self.dx/2,self.xmax-self.dx/2,self.nx)
        self.y=linspace(self.ymin+self.dy/2,self.ymax-self.dy/2,self.ny)
        self.zp=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)
        self.z=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)        


        def ZP(z):
            n = 10.0
            return - 100.0 * (n - 1.0) * pow(z / 100.0, (n + 1.0)/n) + n * z

        def DZPDz(z):
            n = 10.0
            return - (n - 1.0) * (n + 1.0) / n * pow(z / 100.0, 1.0 / n) + n

        def Z(zp):
            # Inverts ZP(z) using Newton's method
            i = 0
            error = 1.0
            lastroot = zp
            while (error>1e-12):
                newroot = lastroot - (ZP(lastroot) - zp) / DZPDz(lastroot)
                error = fabs((newroot - lastroot) / lastroot)
                lastroot=newroot
                i=i+1
                if (i > 100):
                    print "Iterations exceeded in Z(zprime). Quitting."
                    return newroot
            return newroot

        for k in range(self.nz):
            self.z[k] = Z(self.zp[k])

        self.phi=array(hdfphi[hdfphi.items()[0][0]])
        hdfrho = h5py.File(rhofile,'r')
        self.rho=array(hdfrho[hdfrho.items()[0][0]])
        if LogEField == 1:
            Exfile = dir+'/'+filebase+'_'+str(run)+'_Ex'
            Eyfile = dir+'/'+filebase+'_'+str(run)+'_Ey'
            Ezfile = dir+'/'+filebase+'_'+str(run)+'_Ez'
            hdfEx = h5py.File(Exfile,'r')
            self.Ex=array(hdfEx[hdfEx.items()[0][0]])
            hdfEy = h5py.File(Eyfile,'r')
            self.Ey=array(hdfEy[hdfEy.items()[0][0]])
            hdfEz = h5py.File(Ezfile,'r')
            self.Ez=array(hdfEz[hdfEz.items()[0][0]])

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

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = sys.argv[2]
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]

dat = Array3dHDF5(outputfiledir, outputfilebase, ConfigData["LogEField"], run) 

# This holds all of the data
ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixel = ConfigData["GridsPerPixel"]
nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

# A couple of things to customize the plots
PlotEField = bool(ConfigData["PlotEField"])
EdgePlot = bool(ConfigData["EdgePlot"])

nxcenter = nxx/2
nycenter = nyy/2
nxcenter2 = nxcenter

if EdgePlot:
    NumPixelsPlotted = 8 # Special for edge plots
    nycenter2 = nycenter - 16 * GridsPerPixel * ScaleFactor
    nymax = nycenter + 21 * ScaleFactor * GridsPerPixel
    nymin = nycenter - 21 * ScaleFactor * GridsPerPixel

else:
    NumPixelsPlotted = 4
    nycenter2 = nycenter
    nymin = nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2
    nymax = nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2

nxmin = nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2
nxmax = nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2

nzmin = 0
nzmax = 16 * ScaleFactor

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})

print "Making array edge potential plots\n"
figure()
suptitle("Array Edge Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplot(2,2,1)
title("Front Edge")
ylim(-20.0, 10.0)
for slicez in [0,1,2,3,10]:
    plot(dat.x[:],dat.phi[:,0,slicez], label = 'z=%.1f'%dat.z[slicez])
legend()

subplot(2,2,2)
title("Back Edge")
ylim(-20.0, 10.0)
for slicez in [0,1,2,3,10]:
    plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = 'z=%.1f'%dat.z[slicez])
legend()

subplot(2,2,3)
title("Left Edge")
if EdgePlot:
    ylim(-75.0, 25.0)
else:
    ylim(-20.0, 10.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[0,:,slicez], label = 'z=%.1f'%dat.z[slicez])
legend()

subplot(2,2,4)
title("Right Edge")
if EdgePlot:
    ylim(-75.0, 25.0)
else:
    ylim(-20.0, 10.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = 'z=%.1f'%dat.z[slicez])
legend()

savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials.pdf")

print "Making 1D potential and Charge Density plots\n"
figure()

suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 100
numzs = 100


subplot(2,3,1)
title("Phi-Collect Gate")
plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0,(dat.y[nycenter2]+dat.y[nycenter2-1])/2.0))

nxcenter3 = nxcenter2 + 4 * GridsPerPixel*ScaleFactor

plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter2]+dat.y[nycenter2-1])/2.0))
legend(loc = "lower left")
xlabel("Z-Dimension (microns)")
ylim(-10.0, 15.0)
xlim(0.0,4.0)

subplot(2,3,4)
title("Rho-Collect Gate")
zs = []
rhos = []
for i in range(1,numzs):
    for m in [-1,0,1]:
        zs.append(dat.z[i] - 0.4 * (dat.z[i] - dat.z[i+m]))
        rhos.append((dat.rho[nxcenter2,nycenter2,i]+dat.rho[nxcenter2-1,nycenter2,i]+dat.rho[nxcenter2,nycenter2-1,i]+dat.rho[nxcenter2-1,nycenter2-1,i])/4.0)

plot(zs, rhos, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0,(dat.y[nycenter2]+dat.y[nycenter2-1])/2.0))

zs = []
rhos = []
for i in range(1,numzs):
    for m in [-1,0,1]:
        zs.append(dat.z[i] - 0.4 * (dat.z[i] - dat.z[i+m]))
        rhos.append((dat.rho[nxcenter3,nycenter2,i]+dat.rho[nxcenter3-1,nycenter2,i]+dat.rho[nxcenter3,nycenter2-1,i]+dat.rho[nxcenter3-1,nycenter2-1,i])/4.0)

plot(zs, rhos, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter2]+dat.y[nycenter2-1])/2.0))
legend(loc = "lower left")
xlabel("Z-Dimension (microns)")
ylim(-60.0, 40.0)
xlim(0.0,4.0)
nxcenter3 = nxcenter2 + 3 * GridsPerPixel * ScaleFactor / 2
nycenter3 = nycenter2
nycenter4 = nycenter2 + GridsPerPixel * ScaleFactor / 2
subplot(2,3,2)
title("Phi-ChanStop")
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter3]+dat.y[nycenter3-1])/2.0))
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter4,0:phinumzs]+dat.phi[nxcenter3-1,nycenter4,0:phinumzs]+dat.phi[nxcenter3,nycenter4-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter4-1,0:phinumzs])/4.0, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter4]+dat.y[nycenter4-1])/2.0))
legend(loc = "lower left")
xlabel("Z-Dimension (microns)")
ylim(-20.0, 15.0)
xlim(0.0,10.0)
subplot(2,3,5)
title("Rho-ChanStop")
zs = []
rhos = []
for i in range(1,numzs):
    for m in [-1,0,1]:
        zs.append(dat.z[i] - 0.4 * (dat.z[i] - dat.z[i+m]))
        rhos.append((dat.rho[nxcenter3,nycenter3,i]+dat.rho[nxcenter3-1,nycenter3,i]+dat.rho[nxcenter3,nycenter3-1,i]+dat.rho[nxcenter3-1,nycenter3-1,i])/4.0)
plot(zs, rhos, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter3]+dat.y[nycenter3-1])/2.0))

#print "1"
#print dat.rho[nxcenter3, nycenter3,0:phinumzs]
#print dat.phi[nxcenter3, nycenter3,0:phinumzs]
#print dat.rho[nxcenter3-GridsPerPixel:nxcenter3+GridsPerPixel,nycenter3,5]
#print dat.rho[nxcenter3-GridsPerPixel:nxcenter3+GridsPerPixel,nycenter3,2]
#print dat.phi[nxcenter3-GridsPerPixel:nxcenter3+GridsPerPixel,nycenter3,5]

zs = []
rhos = []
for i in range(1,numzs):
    for m in [-1,0,1]:
        zs.append(dat.z[i] - 0.4 * (dat.z[i] - dat.z[i+m]))
        rhos.append((dat.rho[nxcenter3,nycenter4,i]+dat.rho[nxcenter3-1,nycenter4,i]+dat.rho[nxcenter3,nycenter4-1,i]+dat.rho[nxcenter3-1,nycenter4-1,i])/4.0)

plot(zs, rhos, label = "x = %.2f, y = %.2f"%((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0,(dat.y[nycenter4]+dat.y[nycenter4-1])/2.0))

#print "2"
#print dat.rho[nxcenter3, nycenter4,0:phinumzs]
#print dat.phi[nxcenter3, nycenter3,0:phinumzs]


legend(loc = "lower left")
xlabel("Z-Dimension (microns)")
ylim(-200.0, 200.0)
xlim(0.0,10.0)

nxcenter3 = nxcenter2
nycenter3 = nycenter2 + GridsPerPixel * ScaleFactor / 2
subplot(2,3,3)
title("Phi-Barrier Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0),((dat.y[nycenter3]+dat.y[nycenter3-1])/2.0)))
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-20.0, 15.0)
xlim(0.0,4.0)
subplot(2,3,6)
title("Rho-Barrier Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter3]+dat.x[nxcenter3-1])/2.0),((dat.y[nycenter3]+dat.y[nycenter3-1])/2.0)))
zs = []
rhos = []
for i in range(1,numzs):
    for m in [-1,0,1]:
        zs.append(dat.z[i] - 0.4 * (dat.z[i] - dat.z[i+m]))
        rhos.append((dat.rho[nxcenter3,nycenter3,i]+dat.rho[nxcenter3-1,nycenter3,i]+dat.rho[nxcenter3,nycenter3-1,i]+dat.rho[nxcenter3-1,nycenter3-1,i])/4.0)

plot(zs, rhos)
xlabel("Z-Dimension (microns)")
ylim(-40.0, 40.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials.pdf")

print "Making 1D potential Plots #2 \n"
figure()

suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
subplot(1,2,1, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-20.0, 20.0, 21)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
#colorbar()

subplot(1,2,2)
title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
ylim(-10.0, 20.0)
xlim(dat.x[nxmin],dat.x[nxmax])
xlabel("X,Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_2.pdf")

print "Making 1D potential Plots #3 \n"
figure()

suptitle("1D Potentials in Isolation Regions. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
nxcenter3 = nxcenter2 + GridsPerPixel * ScaleFactor / 2
nycenter3 = nycenter2 + GridsPerPixel * ScaleFactor / 2
subplot(1,2,1, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-20.0, 20.0, 21)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter3],dat.y[nycenter3]],ls = "-", color="k")
plot([dat.x[nxcenter3],dat.x[nxcenter3]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
#colorbar()

subplot(1,2,2)
title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter3, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter3])
plot(dat.y[nymin:nymax],dat.phi[nxcenter3,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter3])
ylim(-10.0, 10.0)
xlim(dat.x[nxmin],dat.x[nxmax])
xlabel("X,Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_3.pdf")


print "Making summary plots\n"
figure()
suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(2,2,1, aspect = 1)
title("Phi, z = 0.0")
if EdgePlot:
    levels = linspace(-40.0, 10.0, 51)
else:
    levels = linspace(-10.0, 10.0, 21)
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
#imshow(dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,interpolation = 'None')
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

subplot(2,2,2, aspect = 1)
nxmin2 = nxmin - 9 * GridsPerPixel * ScaleFactor / 2
nymin2 = nymin - 9 * GridsPerPixel * ScaleFactor / 2
rho0 = dat.rho[nxmin2,nymin2,slicez+1]

title("Rho, z = %.2f - %.2f"%(dat.z[nzmin],dat.z[nzmax]))
levels = linspace(-10.0,10.0,41)
plotarray = array(dat.rho[nxmin:nxmax,nymin:nymax,nzmin:nzmax].sum(axis=2)/(nzmax-nzmin))
contour(xx,yy,plotarray, levels, lw=0.1)
contourf(xx,yy,plotarray, levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

slicez = 8 * ScaleFactor
subplot(2,2,3, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
if EdgePlot:
    levels = linspace(-40.0, 20.0, 61)
else:
    levels = linspace(-20.0, 20.0, 21)
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
colorbar()

slicez = 16 * ScaleFactor
subplot(2,2,4, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1.pdf")

figure()
suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)

if EdgePlot:
    levels = linspace(-40.0, 10.0, 51)
else:
    levels = linspace(-20.0, 10.0, 31)

subplot(1,2,1)
title("Phi and (-)E in Gate Region. y = %.2f"%dat.y[nycenter2])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 16 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels,lw=0.1)
contourf(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels)
colorbar()
if ConfigData["LogEField"] == 1 and PlotEField:
    nzmin = 1
    [zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
    quiver(xx, zz, dat.Ex[nxmin:nxmax,nycenter2,nzmin:nzmax], dat.Ez[nxmin:nxmax,nycenter2,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])


subplot(1,2,2)
title("Phi and (-)E in Gate Region. x = %.2f"%dat.x[nxcenter2])
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contour(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels,lw=0.1)
contourf(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels)
colorbar()
if ConfigData["LogEField"] == 1 and PlotEField:
    nzmin = 1
    [zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
    quiver(yy, zz, dat.Ey[nxcenter2,nymin:nymax,nzmin:nzmax], dat.Ez[nxcenter,nymin:nymax,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2.pdf")

if ConfigData["LogPixels"] == 1:
    # Next, plots of the pixel boundaries
    print "Making pixel plots\n"
    figure()
    rcParams['contour.negative_linestyle'] = 'solid'
    #rcParams.update({'font.size': 18})

    Channelkmax = int(ConfigData["ChannelDepth"] / (ConfigData["PixelSize"] / ConfigData["GridsPerPixel"]))
    suptitle("CCD Pixel Plots. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 24)
    plotcounter = 1
    subplots_adjust(hspace=0.3, wspace=0.1)

    filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_Pts"
    file = open(filename,"r")
    lines = file.readlines()
    file.close()
    if len(lines) < 2:
        print "No data in Pts file.  Quitting"
        sys.exit()
    redsx=[]
    redsy=[]
    blacksx=[]
    blacksy=[]
    plottedxin = -1000.0
    plottedyin = -1000.0
    lines.remove(lines[0])
    for line in lines:
        values = line.split()
        zout = float(values[5])    
        if zout > dat.z[Channelkmax + 3]:
            continue # Skip these in case of LogTracePaths = 1
        xin = float(values[0])
        yin = float(values[1])
        xout = float(values[3])
        yout = float(values[4])
        if isnan(xout) or isnan(yout):
            print "xin = %.3f, yin = %.3f is a nan"
            continue
        pixxout = int(xout/10.0)
        pixyout = int(yout/10.0)

        if xin>19.5 and xin < 20.5 and yin >19.5 and yin<20.5:
            print xin, yin, pixxout, pixyout, plottedxin, plottedyin

        if (xin > plottedxin - .001) and (xin < plottedxin + .001) and \
           (yin > plottedyin - .001) and (yin < plottedyin + .001):
            continue # Skips it if it is already plotted
        if (pixxout + pixyout) % 2 == 0:
            redsx.append(xin)
            redsy.append(yin)
        else:
            blacksx.append(xin)
            blacksy.append(yin)

        plottedxin = xin
        plottedyin = yin

    subplot(1,1,1,aspect=1)
    title("Pixel Boundaries",fontsize = 12)
    if ConfigData["PixelBoundaryTestType"] == 0:
        spotsize = 10.0 * ConfigData["PixelBoundaryStepSize"][0] * ConfigData["PixelBoundaryStepSize"][1]
    else:
        spotsize = 0.1
    scatter(redsx,redsy,s=spotsize,color="red")
    scatter(blacksx,blacksy,s=spotsize,color="black")

    if EdgePlot:
        for linex in linspace(120.0,200.0,9):
            plot((linex,linex),(20.0,70.0),linewidth=1.0, color='blue')

    xlabel("X(microns)",fontsize = 18)
    ylabel("Y(microns)",fontsize = 18)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
    ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


    savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels.pdf")

if ConfigData["LogPixelPaths"] == 1 and ConfigData["LogPixels"] == 1:
    # Last, plots of the electron paths
    print "Making array electron path plots\n"

    yline = (ConfigData["PixelBoundaryLowerLeft"][1] + ConfigData["PixelBoundaryUpperRight"][1] + ConfigData["PixelBoundaryStepSize"][1]) / 2.0
    xline = (ConfigData["PixelBoundaryLowerLeft"][0] + ConfigData["PixelBoundaryUpperRight"][0] + ConfigData["PixelBoundaryStepSize"][0]) / 2.0

    figure()
    vertical_zoom = 1
    suptitle("Electron Path Plot - Vertical Zoom = %d"%vertical_zoom, fontsize = 24)
    subplots_adjust(wspace=0.2)
    subplot(1,2,1,aspect=vertical_zoom)
    oldxin = 1000000.0
    lines.remove(lines[0])
    for line in lines:
        values = line.split()
        xin = float(values[0])
        yin = float(values[1])
        xout = float(values[3])
        yout = float(values[4])
        zout = float(values[5])
        if (yin < yline - .10) or (yin > yline + .10):
            continue


        if abs(xin - oldxin) > .001:
            if oldxin < 100000.0:
                pixxin = int(oldxin/10.0)
                if pixxin % 2 == 0:
                    color = "red"
                else:
                    color = "black"
                plot(xpaths, zxpaths, color = color, linewidth = 0.1)
            oldxin = xin
            xpaths=[]
            zxpaths=[]
            
        if isnan(xout) or isnan(yout) or isnan(zout):
            continue
        xpaths.append(xout)
        zxpaths.append(zout)

    ylabel("Z(microns)")
    xlabel("X (microns)")
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])

    subplot(1,2,2,aspect=vertical_zoom)
    oldyin = 1000000.0
    for line in lines:
        values = line.split()
        xin = float(values[0])
        yin = float(values[1])
        xout = float(values[3])
        yout = float(values[4])
        zout = float(values[5])

        if (xin < xline - .10) or (xin > xline + .10):
            continue

        if abs(yin - oldyin) > .001:
            if oldyin < 100000.0:
                pixyin = int(oldyin/10.0)
                if pixyin % 2 == 0:
                    color = "red"
                else:
                    color = "black"

                plot(ypaths, zypaths, color = color, linewidth = 0.1)
            oldyin = yin
            ypaths=[]
            zypaths=[]
            
        if isnan(xout) or isnan(yout) or isnan(zout):
            continue
        ypaths.append(yout)
        zypaths.append(zout)

    ylabel("Z(microns)")
    xlabel("Y (microns)")
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])
    savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths.pdf")
