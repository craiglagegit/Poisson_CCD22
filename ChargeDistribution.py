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
    def __init__(self, dir, filebase, LogEField):
        phifile = dir+'/'+filebase+'_phi'
        rhofile = dir+'/'+filebase+'_rho'
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
        self.z=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)

        self.phi=array(hdfphi[hdfphi.items()[0][0]])
        hdfrho = h5py.File(rhofile,'r')
        self.rho=array(hdfrho[hdfrho.items()[0][0]])
        if LogEField == 1:
            Exfile = dir+'/'+filebase+'_Ex'
            Eyfile = dir+'/'+filebase+'_Ey'
            Ezfile = dir+'/'+filebase+'_Ez'
            hdfEx = h5py.File(Exfile,'r')
            self.Ex=array(hdfEx[hdfEx.items()[0][0]])
            hdfEy = h5py.File(Eyfile,'r')
            self.Ey=array(hdfEy[hdfEy.items()[0][0]])
            hdfEz = h5py.File(Ezfile,'r')
            self.Ez=array(hdfEz[hdfEz.items()[0][0]])

class Array3dHDF5Elec(object):
    def __init__(self, dir, filebase, n):
        elecfile = dir+'/'+filebase+'_'+str(n)+'_Elec'
        hdfelec = h5py.File(elecfile,'r')
        Dimension = hdfelec[hdfelec.items()[0][0]].attrs[u'Dimension']
        self.nx=Dimension[0]
        self.ny=Dimension[1]
        self.nz=Dimension[2]
        
        Lower_Left = hdfelec[hdfelec.items()[0][0]].attrs[u'Lower_Left']
        self.xmin=Lower_Left[0]
        self.ymin=Lower_Left[1]
        self.zmin=Lower_Left[2]

        Upper_Right = hdfelec[hdfelec.items()[0][0]].attrs[u'Upper_Right']
        self.xmax=Upper_Right[0]
        self.ymax=Upper_Right[1]
        self.zmax=Upper_Right[2]
        
        self.dx=(self.xmax-self.xmin)/self.nx
        self.dy=(self.ymax-self.ymin)/self.ny
        self.dz=(self.zmax-self.zmin)/self.nz
        self.volume = self.dx * self.dy * self.dz
        
        self.x=linspace(self.xmin+self.dx/2,self.xmax-self.dx/2,self.nx)
        self.y=linspace(self.ymin+self.dy/2,self.ymax-self.dy/2,self.ny)
        self.z=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)

        self.elec=array(hdfelec[hdfelec.items()[0][0]])

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
run = int(sys.argv[2])

ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]


# This holds all of the data
ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixel = ConfigData["GridsPerPixel"]
dat = Array3dHDF5Elec(outputfiledir, outputfilebase, run) 
nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1
NumPixelsPlotted = 3
nxcenter = nxx/2
nycenter = nyy/2
nxmin = nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2
nxmax = nxmin + (NumPixelsPlotted * ScaleFactor * GridsPerPixel)
nymin = nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixel)/2
nymax = nymin + (NumPixelsPlotted * ScaleFactor * GridsPerPixel)
levels = linspace(-3.0, 7.0, 101)

Total_elec = dat.elec.sum()

print "Total Electrons = ",Total_elec

fig = figure(figsize = (12,12))
suptitle("Well Filling with Two Collecting Phases", fontsize = 36)

ax1=axes([0.10,0.40,0.50,0.50],aspect=1)
ax1.set_title("X-Y Slice")
ax1.set_xticks([])
ax1.set_yticks([])
plotarray = sum(dat.elec[nxmin:nxmax,nymin:nymax,:],axis = 2)+0.1
for nx in range(ScaleFactor*GridsPerPixel,nxmax-nxmin,ScaleFactor*GridsPerPixel):
    for ny in range(0,nymax-nymin):
        plotarray[nx,ny] = 0.001
for ny in range(ScaleFactor*GridsPerPixel,nymax-nymin,ScaleFactor*GridsPerPixel):
    for nx in range(0,nxmax-nxmin):
        plotarray[nx,ny] = 0.001
for nx in range(ScaleFactor*GridsPerPixel-1,nxmax-nxmin,ScaleFactor*GridsPerPixel):
    for ny in range(0,nymax-nymin):
        plotarray[nx,ny] = 0.001
for ny in range(ScaleFactor*GridsPerPixel-1,nymax-nymin,ScaleFactor*GridsPerPixel):
    for nx in range(0,nxmax-nxmin):
        plotarray[nx,ny] = 0.001

ax1.imshow(log10(transpose(plotarray)), interpolation = 'nearest')

ax2=axes([0.10,0.10,0.50,0.20])
ax2.set_title("X-Z Slice")
ax2.set_xticks([])
ax2.set_yticks([])
plotarray = sum(dat.elec[nxmin:nxmax,:,0:20],axis = 1)+0.1
for nx in range(0,nxmax-nxmin):
    plotarray[nx,0] = 0.001

ax2.imshow(log10(transpose(fliplr(plotarray))), interpolation = 'nearest')

ax3=axes([0.70,0.40,0.20,0.50])
ax3.set_title("Y-Z Slice")
ax3.set_xticks([])
ax3.set_yticks([])
plotarray = sum(dat.elec[:,nymin:nymax,0:20],axis = 0)+0.1
for ny in range(0,nymax-nymin):
    plotarray[ny,0] = 0.001

ax3.imshow(log10(fliplr(plotarray)), interpolation = 'nearest')
savefig(outputfiledir+"/plots/ChargeDistribution_%d.pdf"%run)
close(fig)






