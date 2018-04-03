#!/usr/bin/env python

#Author: Craig Lage, UCDavis
#Date: 2-Apr-18

#This program manages the running of multiple Poisson tests
from pylab import *
import os, sys, time, subprocess
#from mpi4py import MPI
import mpi
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *
#****************MAIN PROGRAM*****************
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()

size = mpi.size
rank = mpi.rank
run = rank

NumRuns = size

newdir = subprocess.Popen('mkdir -p data/run_%d'%run, shell=True) 
subprocess.Popen.wait(newdir)
incfgfile = 'data/run/pixel.cfg'
outcfgfile = 'data/run_%d/pixel.cfg'%run
ConfigData = ReadConfigFile(incfgfile)
New_Fe55_Cfg_File(incfgfile, outcfgfile, run)
time.sleep(1.0)

cmd = '~/Software/Poisson_CCD_Hole20/src_fe55_2apr18/Poisson data/run_%d/pixel.cfg'%run
print "Launching %s"%cmd
cpp_job = subprocess.Popen(cmd, shell=True)
subprocess.Popen.wait(cpp_job)
time.sleep(1.0)



