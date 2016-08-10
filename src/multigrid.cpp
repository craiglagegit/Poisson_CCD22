/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** multigrid.cpp **************************

#include "multigrid.h"

MultiGrid::MultiGrid(string inname) //Constructor                                                                                            
{
  // This reads in the data from the poisson.cfg
  // file, sets the initial conditions, and solves
  // Poisson's equation using SOR and Multi-Grid methods
  double setup_time, solution_time, efield_time, trace_time;
  time_t time1, time2;
  time1 = time(NULL);
  
  //Set the random number seed
  srand48( (unsigned int) time(NULL));

  // First we read in the configuration information
  ReadConfigurationFile(inname);
  printf("Finished Reading config file\n");
  // Then, we build the multigrid arrays and set the initial conditions
  nsteps = 4 + (int)(log2(ScaleFactor));
  // nsteps is the number of reduction steps in the VCycle.
  // This gives 4 grid cells in the z-direction at the coarsest scale
  phi = new Array3D*[nsteps+1];
  rho = new Array3D*[nsteps+1];
  elec = new Array3D*[1];
  hole = new Array3D*[1];    
  BCType = new Array2D*[nsteps+1];  
  E = new Array3D*[3];

  BuildArrays(phi, rho, elec, E, BCType);
  printf("Finished Building Arrays. \n");
  Channelkmin = phi[0]->ZIndex(GateOxide * EPSILON_SI / EPSILON_OX) + 1;
  ChannelStopkmin = Channelkmin;

  SetInitialVoltages(phi[0], BCType[0]);
  SetFixedCharges(rho[0], BCType[0]); // Place fixed charges      
  SetInitialHoles(rho[0], hole[0]);
  SetInitialElectrons(rho[0], elec[0]);  // This sets initial electrons
  time2 = time(NULL);
  setup_time = difftime(time2, time1);
  printf("Finished Setting ICs. Setup time = %.3f seconds\n",setup_time);
  fflush(stdout);

  // Now we run NumSteps cycle, adding NumElec electrons each step and re-solving
  // Poisson's equation at each step.
  int m, kmax;
  string StepNum;
  string underscore = "_";
  for (m=0; m<NumSteps; m++)
    {
      time1 = time(NULL);
      StepNum = boost::lexical_cast<std::string>(m);      
      if (m == 0) kmax = 20;
      else kmax = 1;

      for (int k=0; k<kmax; k++)
	{
	  SetFixedCharges(rho[0], BCType[0]); // Place fixed charges    
	  if (!(m == 0 && k == 0))
	    {
	      AdjustHoles(phi[0], rho[0], hole[0]);
	    }
	  FillRho(rho[0], elec[0], hole[0]); // Add mobile charges.
	  // Now we cycle through the VCycles to solve Poisson's equation
	  for (int n=0; n<iterations; n++)
	    {
	      if (k == kmax - 1)
		{
		  VCycle(phi, rho, BCType, w, nsteps, ncycle);
		}
	      else
		{
		  VCycle(phi, rho, BCType, w, nsteps, 4);
		}
	    }
	  time2 = time(NULL);
	  solution_time = difftime(time2, time1);
	  printf("Finished solving Poisson's equation. Solution time = %.3f seconds\n",solution_time);
	  time1 = time(NULL);
	  Gradient(phi[0], E);
	}
      // Next we calculate the E fields
      Gradient(phi[0], E);
      time2 = time(NULL);
      efield_time = difftime(time2, time1);
      printf("Finished calculating E Fields. E Field time = %.3f seconds\n",efield_time);
      time1 = time(NULL);

      // Now we trace the electrons.
      if (PixelAreas >= 0 && m % PixelAreas == 0)
	{
	  CalculatePixelAreas(m);
	}
      if (PixelBoundaryTestType == 0)
	{
	  TraceGrid(m);
	}
      if (PixelBoundaryTestType == 1)
	{
	  TraceSpot(m);
	}
      if (PixelBoundaryTestType == 2)
	{
	  TraceRegion(m);
	}
      if (PixelBoundaryTestType == 3)
	{
	  TraceMultipleSpots(m);
	}
      // Now, we write out the potential and charge density results

      time2 = time(NULL);
      trace_time = difftime(time2, time1);
      printf("Finished tracing electrons. Trace time = %.3f seconds\n",trace_time);
      
      if (SaveElec !=0 && m % SaveElec == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);	  
	}

      if (SaveData !=0 && m % SaveData == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "phi", phi[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "rho", rho[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ex", E[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ey", E[1]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ez", E[2]);
	}
    }
  return;
}

MultiGrid::~MultiGrid() //Destructor                                                                                            
{
  int n;
  for (n=0; n<nsteps+1; n++)
    {
      delete phi[n];
      delete rho[n];
      delete BCType[n];      
    }
  for (n=0; n<3; n++)
    {
      delete E[n];
    }
  delete[] E;
  delete[] phi;
  delete[] rho;
  delete elec[0];
  delete[] elec;
  delete[] BCType;
  return;
}

void MultiGrid::ReadConfigurationFile(string inname)
{
  // Poisson solver constants
  w = GetDoubleParam(inname, "w", 1.9);			// Successive Over-Relaxation factor
  ncycle = GetIntParam(inname, "ncycle", 100);		// Number of SOR sysles at each resolution
  iterations =  GetIntParam(inname, "iterations", 3);	// Number of VCycles

  // Overall setup
  
  SaveData =  GetIntParam(inname, "SaveData", 1);     // 0 - Save only Pts, N save phi,rho,E every Nth step
  SaveElec =  GetIntParam(inname, "SaveElec", 1);     // 0 - Save only Pts, N save Elec every Nth step  
  ScaleFactor =  GetIntParam(inname, "ScaleFactor", 1);     // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 5/6 micron, 128 grids in the z-direction
  PixelSize = GetDoubleParam(inname, "PixelSize", 10.0);    // Pixel size in microns
  GridsPerPixel = GetIntParam(inname, "GridsPerPixel", 12); // Grids per pixel at ScaleFactor = 1
  GridsPerPixel = GridsPerPixel * ScaleFactor;
  GridSpacing = PixelSize / (double)GridsPerPixel;

  Nx = GetIntParam(inname, "Nx", 160);                // Number of grids in x at ScaleFactor = 1
  Nx = Nx * ScaleFactor;
  Ny = GetIntParam(inname, "Ny", 160);                // Number of grids in y at ScaleFactor = 1
  Ny = Ny * ScaleFactor;
  Nz = 160;                                           // Number of grids in z at ScaleFactor = 1
  Nz = Nz * ScaleFactor;
  Nzelec = 32;                                        // Number of grids in z in electron array at ScaleFactor = 1
  Nzelec = Nzelec * ScaleFactor;
  XBCType = GetIntParam(inname, "XBCType", 1);        // 0 - Free BC, 1 - Periodic BC
  YBCType = GetIntParam(inname, "YBCType", 1);        // 0 - Free BC, 1 - Periodic BC  
  SimulationRegionLowerLeft = new double[2];
  SimulationRegionLowerLeft[0] = 0.0;  SimulationRegionLowerLeft[1] = 0.0;
  SimulationRegionLowerLeft = GetDoubleList(inname, "SimulationRegionLowerLeft", 2, SimulationRegionLowerLeft);  
  // Voltages and Charges
  Vbb = GetDoubleParam(inname, "Vbb", -50.0);		        // Back bias
  Vparallel_lo = GetDoubleParam(inname, "Vparallel_lo", -8.0);	// Parallel Low Voltage
  Vparallel_hi = GetDoubleParam(inname, "Vparallel_hi", 4.0);	// Parallel High Voltage
  Vserial_lo = GetDoubleParam(inname, "Vserial_lo", -6.0);	// Serial Low Voltage
  Vaverage = (8.0 * Vparallel_lo + 4.0 * Vparallel_hi) / 12.0;
  Vchannelstop = GetDoubleParam(inname, "Vchannelstop", 0.0);
  Vscupper = GetDoubleParam(inname, "Vscupper", 5.0);
  GateOxide = GetDoubleParam(inname, "GateOxide", 0.15);
  BackgroundDoping = GetDoubleParam(inname, "BackgroundDoping", -1.0E12);
  ChannelStopDoping = GetDoubleParam(inname, "ChannelStopDoping", -1.0E12);
  ChannelStopDepth = GetDoubleParam(inname, "ChannelStopDepth", 1.0);
  ChannelStopProfile = GetIntParam(inname, "ChannelStopProfile", 0);  
  ChannelStopWidth = GetDoubleParam(inname, "ChannelStopWidth", 1.0);
  ChannelDoping = GetDoubleParam(inname, "ChannelDoping", -5.0E11);
  ChannelDepth = GetDoubleParam(inname, "ChannelDepth", 1.0);
  ChannelProfile = GetIntParam(inname, "ChannelProfile", 0);    
  UndepletedChannelStop = GetIntParam(inname, "UndepletedChannelStop", 0);

  // Pixel Regions
  int i, j, k;
  string regionnum, fillednum;
  NumberofPixelRegions = GetIntParam(inname, "NumberofPixelRegions", 0);
  PixelRegionLowerLeft = new double*[NumberofPixelRegions];
  PixelRegionUpperRight = new double*[NumberofPixelRegions];
  NumberofFilledWells = new int[NumberofPixelRegions];
  CollectingPhases = GetIntParam(inname, "CollectingPhases", 1);
  FilledPixelCoords = new double**[NumberofPixelRegions];
  CollectedCharge = new int*[NumberofPixelRegions];
  CollectedChargeZmin = GetDoubleParam(inname, "CollectedChargeZmin", 0.60);
  CollectedChargeZmax = GetDoubleParam(inname, "CollectedChargeZmax", 0.90);  
  CollectedChargeXmin = GetDoubleParam(inname, "CollectedChargeXmin", 1.0);
  CollectedChargeXmax = GetDoubleParam(inname, "CollectedChargeXmax", 9.0);  
  CollectedChargeYmin = GetDoubleParam(inname, "CollectedChargeYmin", 1.66);
  CollectedChargeYmax = GetDoubleParam(inname, "CollectedChargeYmax", 8.34);  

  for (i=0; i<NumberofPixelRegions; i++)
    {
      PixelRegionLowerLeft[i] = new double[2];
      PixelRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  PixelRegionLowerLeft[i][j] = 0.0;
	  PixelRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofPixelRegions; i++)
    {
      regionnum = boost::lexical_cast<std::string>(i);
      PixelRegionLowerLeft[i] = GetDoubleList(inname, "PixelRegionLowerLeft_"+regionnum, 2, PixelRegionLowerLeft[i]);
      PixelRegionUpperRight[i] = GetDoubleList(inname, "PixelRegionUpperRight_"+regionnum, 2, PixelRegionUpperRight[i]);
      NumberofFilledWells[i] = GetIntParam(inname, "NumberofFilledWells_"+regionnum, 0);
      CollectedCharge[i] = new int[NumberofFilledWells[i]];
      FilledPixelCoords[i] = new double*[NumberofFilledWells[i]];
      for (j=0; j<NumberofFilledWells[i]; j++)
	{
	  fillednum = boost::lexical_cast<std::string>(j);
	  CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0);
	  FilledPixelCoords[i][j] = new double[2];
	  for (k=0; k<2; k++)
	    {
	      FilledPixelCoords[i][j][k] = 0.0;
	    }
	  FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j]);
	}
    }
  
  // Fixed Voltage Regions
  NumberofFixedRegions = GetIntParam(inname, "NumberofFixedRegions", 0);
  FixedRegionLowerLeft = new double*[NumberofFixedRegions];
  FixedRegionUpperRight = new double*[NumberofFixedRegions];
  FixedRegionVoltage = new double[NumberofFixedRegions];
  FixedRegionDoping = new double[NumberofFixedRegions];
  FixedRegionChargeDepth = new double[NumberofFixedRegions];
  FixedRegionBCType = new int[NumberofFixedRegions];  

  for (i=0; i<NumberofFixedRegions; i++)
    {
      FixedRegionLowerLeft[i] = new double[2];
      FixedRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  FixedRegionLowerLeft[i][j] = 0.0;
	  FixedRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofFixedRegions; i++)
    {
      regionnum = boost::lexical_cast<std::string>(i);
      FixedRegionLowerLeft[i] = GetDoubleList(inname, "FixedRegionLowerLeft_"+regionnum, 2, FixedRegionLowerLeft[i]);
      FixedRegionUpperRight[i] = GetDoubleList(inname, "FixedRegionUpperRight_"+regionnum, 2, FixedRegionUpperRight[i]);
      FixedRegionVoltage[i] = GetDoubleParam(inname, "FixedRegionVoltage_"+regionnum,0.0);
      FixedRegionDoping[i] = GetDoubleParam(inname, "FixedRegionDoping_"+regionnum,0.0);
      FixedRegionChargeDepth[i] = GetDoubleParam(inname, "FixedRegionChargeDepth_"+regionnum,0.0);
      FixedRegionBCType[i] = GetIntParam(inname, "FixedRegionBCType_"+regionnum,0);      
    }


  // Pixel Boundary Tests

  LogEField = GetIntParam(inname, "LogEField", 0);
  LogPixels = GetIntParam(inname, "LogPixels", 0);
  LogPixelPaths = GetIntParam(inname, "LogPixelPaths", 0);
  PixelAreas = GetIntParam(inname, "PixelAreas", 0);  
  PixelBoundaryTestType = GetIntParam(inname, "PixelBoundaryTestType", 0);    
  PixelBoundaryLowerLeft = new double(2);
  PixelBoundaryUpperRight = new double(2);
  PixelBoundaryStepSize = new double(2);
  for (j=0; j<2; j++)
    {
      PixelBoundaryLowerLeft[j] = 0.0;
      PixelBoundaryUpperRight[j] = 100.0;
      PixelBoundaryStepSize[j] = 1.0;
    }
  PixelBoundaryLowerLeft = GetDoubleList(inname, "PixelBoundaryLowerLeft", 2, PixelBoundaryLowerLeft);
  PixelBoundaryUpperRight = GetDoubleList(inname, "PixelBoundaryUpperRight", 2, PixelBoundaryUpperRight);
  CCDTemperature = GetDoubleParam(inname, "CCDTemperature", 173.0);
  DiffMultiplier = GetDoubleParam(inname, "DiffMultiplier", 1.0);
  SaturationModel = GetIntParam(inname, "SaturationModel", 0);          
  NumDiffSteps = GetIntParam(inname, "NumDiffSteps", 1);        
  NumVertices = GetIntParam(inname,"NumVertices",2);
  ElectronZ0Area = GetDoubleParam(inname,"ElectronZ0Area",100.0);
  ElectronZ0Fill = GetDoubleParam(inname,"ElectronZ0Fill",100.0);

  if (PixelBoundaryTestType == 0)
    {
      NumSteps = GetIntParam(inname, "NumSteps", 100);            
      PixelBoundaryStepSize = new double(2);
      for (j=0; j<2; j++)
	{
	  PixelBoundaryStepSize[j] = 1.0;
	}
      PixelBoundaryStepSize = GetDoubleList(inname, "PixelBoundaryStepSize", 2, PixelBoundaryStepSize);
    }
  if (PixelBoundaryTestType == 1)
    {
      PixelBoundaryNx = GetIntParam(inname, "PixelBoundaryNx", 9);
      PixelBoundaryNy = GetIntParam(inname, "PixelBoundaryNy", 9);  
      NumElec = GetIntParam(inname, "NumElec", 1000);
      NumSteps = GetIntParam(inname, "NumSteps", 100);      
      Sigmax = GetDoubleParam(inname, "Sigmax", 1.0);
      Sigmay = GetDoubleParam(inname, "Sigmay", 1.0);
      Xoffset = GetDoubleParam(inname, "Xoffset", 0.0);
      Yoffset = GetDoubleParam(inname, "Yoffset", 0.0);
    }
  if (PixelBoundaryTestType == 2)
    {
      NumElec = GetIntParam(inname, "NumElec", 1000);
      NumSteps = GetIntParam(inname, "NumSteps", 100);      
    }
  outputfilebase  = GetStringParam(inname,"outputfilebase", "Test"); //Output filename base
  outputfiledir  = GetStringParam(inname,"outputfiledir", "data"); //Output filename directory
  return;
}

void MultiGrid::BuildArrays(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** E, Array2D** BCType)
{
  // Builds the multigrid arrays
  int nx, ny, nz, nxx, nyy, nzz, n;
  double xmin, xmax, ymin, ymax, zmin, zmax, zmaxelec, dx, dy, dz;
  nxx = Nx + 1;
  nyy = Ny + 1;
  nzz = Nz + 1;
  dx = PixelSize / (double)GridsPerPixel;
  dy = dx;
  Xmin = SimulationRegionLowerLeft[0];
  Ymin = SimulationRegionLowerLeft[1];
  Zmin = 0.0;
  Xmax = dx * (double)nxx + Xmin;
  Ymax = dy * (double)nyy + Ymin;
  Zmax = 100.0;
  phi[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz);
  rho[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz);
  zmaxelec = rho[0]->Z(rho[0]->zp[Nzelec] + rho[0]->dzp / 2.0);
  elec[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,Nzelec);
  hole[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,Nzelec);    
  BCType[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);    

  for (n=1; n<nsteps+1; n++)
    {
      nx = (phi[0]->nx - 1) / (int)pow(2,n) + 1;
      ny = (phi[0]->ny - 1) / (int)pow(2,n) + 1;
      nz = (phi[0]->nz - 1) / (int)pow(2,n) + 1;
      dx = phi[0]->dx * (int)pow(2,n);
      dy = phi[0]->dy * (int)pow(2,n);
      dz = phi[0]->dzp * (int)pow(2,n);
      xmin = phi[0]->xmin + phi[0]->dx / 2.0 - dx / 2.0;
      ymin = phi[0]->ymin + phi[0]->dy / 2.0 - dy / 2.0;
      zmin = phi[0]->zmin + phi[0]->dzp / 2.0 - dz / 2.0; 
      xmax = phi[0]->xmax - phi[0]->dx / 2.0 + dx / 2.0;
      ymax = phi[0]->ymax - phi[0]->dy / 2.0 + dy / 2.0;
      zmax = phi[0]->zmax - phi[0]->dzp / 2.0 + dz / 2.0;
      phi[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
      rho[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
      BCType[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);      
    }
  for (n=0; n<3; n++)
    {
      E[n] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz);
    }
  return;
}

void MultiGrid::SetInitialVoltages(Array3D* phi, Array2D* BCType)
{
  int i, j, n, index;
  int PixX, PixY;
  double PixXmin, PixYmin;

  // Potential on top
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  index = i + j * phi->nx + (phi->nz - 1) * phi->nx * phi->ny;
	  phi->data[index] = Vbb;
	}
    }
    // Fixed Potentials or free Boundary Conditions on bottom
  for (n=0; n<NumberofFixedRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  for (j=0; j<phi->ny; j++)
	    {
	      index = i + j * phi->nx;
	      if (phi->x[i] >= FixedRegionLowerLeft[n][0] && phi->x[i] <= FixedRegionUpperRight[n][0] && phi->y[j] >= FixedRegionLowerLeft[n][1] && phi->y[j] <= FixedRegionUpperRight[n][1])
		{
		  if (FixedRegionBCType[n] == 0) //  Fixed potential region
		    {
		      phi->data[index] = FixedRegionVoltage[n];
		    }
		  else // Free Boundary Condition region
		    {
		      BCType->data[index] = 1;
		    }
		}
	    }
	}
    }
    // Potentials in Pixel Region
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  if (phi->x[i] < PixelRegionLowerLeft[n][0] || phi->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<phi->ny; j++)
	    {
	      if (phi->y[j] < PixelRegionLowerLeft[n][1] || phi->y[j] > PixelRegionUpperRight[n][1])
		{ 
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * phi->nx;
	      PixX = (int)((phi->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)((phi->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
	      // Set the gate voltages
	      if (CollectingPhases == 2) // Two collecting phases
		{
		  if (phi->y[j] >= PixYmin + PixelSize/6.0 && phi->y[j] <= PixYmin + 5.0*PixelSize/6.0)
		    {
		      // This is the collection region
		      phi->data[index] = Vparallel_hi;
		    }
		  else
		    {
		      // This is the barrier gate region
		      phi->data[index] = Vparallel_lo;
		    }
		}
	      else if (CollectingPhases == 1) // One collecting phase
		{
		  if (phi->y[j] >= PixYmin + PixelSize/3.0 && phi->y[j] <= PixYmin + 2.0*PixelSize/3.0)
		    {
		      // This is the collection region
		      phi->data[index] = Vparallel_hi;
		    }
		  else
		    {
		      // This is the barrier gate region
		      phi->data[index] = Vparallel_lo;
		    }
		}

	      if ((phi->x[i] <= PixXmin + ChannelStopWidth/2.0 || phi->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0) && UndepletedChannelStop == 1)
		{
		  // This is the Channel Stop Region
		  phi->data[index] = Vchannelstop;
		}
	    }
	}
    }
  printf("Finished setting Boundary Potentials, \n");
  fflush(stdout);
  return;
}

void MultiGrid::SetFixedCharges(Array3D* rho, Array2D* BCType)
{
  int i, j, k, n, index, index2;
  int FixedRegionkmin, FixedRegionkmax;
  FixedRegionkmin = ChannelStopkmin;
  if (ChannelProfile == 0) // Square well
    {
      Channelkmax = rho->ZIndex(rho->Z(rho->zp[Channelkmin]) + ChannelDepth);
    }
  if (ChannelProfile == 1) // Gaussian
    {
      Channelkmax = rho->ZIndex(rho->Z(rho->zp[Channelkmin]) + 4.0 * ChannelDepth);
    }
  if (ChannelStopProfile == 0) // Square well
    {
      ChannelStopkmax = rho->ZIndex(rho->Z(rho->zp[ChannelStopkmin]) + ChannelStopDepth);
    }
  if (ChannelStopProfile == 1) // Gaussian
    {
      ChannelStopkmax = rho->ZIndex(rho->Z(rho->zp[ChannelStopkmin]) + 2.0 * ChannelStopDepth);
    }
  int PixX, PixY;
  double PixXmin, PixYmin, FRChargeDepth, CChargeDepth;
  double ChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
  // ChargeFactor converts doping in cm^-3 into the appropriate units
  double ChannelCharge, FixedRegionCharge;
  ChannelStopCharge = ChannelStopDoping * MICRON_PER_CM  * ChargeFactor;
  ChannelCharge =  ChannelDoping * MICRON_PER_CM  * ChargeFactor;
  CSChargeDepth = rho->Z(rho->zp[ChannelStopkmax] + rho->dzp / 2.0) - rho->Z(rho->zp[ChannelStopkmin] - rho->dzp / 2.0);
  CChargeDepth = rho->Z(rho->zp[Channelkmax] + rho->dzp / 2.0) - rho->Z(rho->zp[Channelkmin] - rho->dzp / 2.0);
  printf("CChargeDepth = %.4f\n",CChargeDepth);
  double ChannelZmin, ChannelZmax, ChannelStopZmin, ChannelStopZmax, ChannelTotal, ChannelStopTotal, DeltaZ;
  ChannelZmin = rho->Z(rho->zp[Channelkmin] - rho->dzp / 2.0);
  ChannelZmax = rho->Z(rho->zp[Channelkmax] + rho->dzp / 2.0);  
  ChannelTotal = erf((ChannelZmax - ChannelZmin) / (sqrt(2.0) * ChannelDepth));
  ChannelStopZmin = rho->Z(rho->zp[ChannelStopkmin] - rho->dzp / 2.0);
  ChannelStopZmax = rho->Z(rho->zp[ChannelStopkmax] + rho->dzp / 2.0);  
  ChannelStopTotal = erf((ChannelStopZmax - ChannelStopZmin) / (sqrt(2.0) * ChannelStopDepth));
  // Set the background charge:

  double Gox_effective = rho->Z(rho->zp[Channelkmin] - rho->dzp / 2.0) * EPSILON_OX / EPSILON_SI;
  printf("In SetFixedCharges, BackgroundCharge = %.4f Effective Gate Oxide thickness = %.3f microns.\n",BackgroundDoping*ChargeFactor, Gox_effective);
  printf("In SetFixedCharges, Channelkmin = %d, Channelkmax = %d, ChannelZWidth = %d grid cells, ChannelCharge = %.4f \n",Channelkmin, Channelkmax, (Channelkmax - Channelkmin + 1), ChannelCharge);
  printf("In SetFixedCharges, ChannelStopkmin = %d, ChannelStopkmax = %d, ChannelStopZWidth = %d grid cells, ChannelStopCharge = %.4f \n",ChannelStopkmin, ChannelStopkmax, (ChannelStopkmax - ChannelStopkmin + 1), ChannelStopCharge);    
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
		{
		  for (k=Channelkmin; k<rho->nz-1; k++)
			{
			  index = i + j * rho->nx + k * rho->nx * rho->ny;
			  rho->data[index] = BackgroundDoping * ChargeFactor;
			}
		}
    }
  
  // Fixed Potentials or free Boundary Conditions on bottom
  for (n=0; n<NumberofFixedRegions; n++)
    {
      rho->ZIndex(rho->Z(rho->zp[Channelkmin]) + ChannelDepth);
      FixedRegionkmax = rho->ZIndex(rho->Z(rho->zp[FixedRegionkmin]) + FixedRegionChargeDepth[n]);      
      FRChargeDepth = rho->Z(rho->zp[FixedRegionkmax] + rho->dzp / 2.0) - rho->Z(rho->zp[FixedRegionkmin] - rho->dzp / 2.0);
      FixedRegionCharge =  FixedRegionDoping[i] * MICRON_PER_CM / FRChargeDepth * ChargeFactor;
      for (i=0; i<rho->nx; i++)
		{
		  for (j=0; j<rho->ny; j++)
			{
			  index = i + j * rho->nx;
			  if (rho->x[i] >= FixedRegionLowerLeft[n][0] && rho->x[i] <= FixedRegionUpperRight[n][0] && rho->y[j] >= FixedRegionLowerLeft[n][1] && rho->y[j] <= FixedRegionUpperRight[n][1])
				{
				  if (FixedRegionBCType[n] == 0) //  Fixed potential region
					{
					  for (k=FixedRegionkmin; k<FixedRegionkmax+1; k++)
						{
						  index2 = index + k * rho->nx * rho->ny;
						  rho->data[index2] = FixedRegionCharge;
						}
					}
				  else // Free Boundary Condition region
					{
					  BCType->data[index] = 1;
					}
				}
			}
		}
    }
  
  // Charges in Pixel Region

  //int nx2 = rho->nx / 2; Debugging
  //int ny2 = rho->ny / 2;  
  //double TotalChannelCharge = 0.0;

  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
		{
		  if (rho->x[i] < PixelRegionLowerLeft[n][0] || rho->x[i] > PixelRegionUpperRight[n][0])
			{
			  continue; // If not in PixelRegion, continue
			}
		  for (j=0; j<rho->ny; j++)
			{
			  if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
				{ 
				  continue; // If not in PixelRegion, continue
				}
			  index = i + j * rho->nx;
			  PixX = (int)((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
			  PixY = (int)((rho->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
			  PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
			  PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
			  // Now set the charges
			  if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0)
				{
				  // This is the Channel Stop Region
				  for (k=Channelkmin; k<ChannelStopkmin; k++)
				    {
				      index2 = index + k * rho->nx * rho->ny;
				      rho->data[index2] = 0.0;
				    }
				  for (k=ChannelStopkmin; k<ChannelStopkmax+1; k++)
					{
					  DeltaZ = rho->Z(rho->zp[k] + rho->dzp / 2.0) - rho->Z(rho->zp[k] - rho->dzp / 2.0);
					  index2 = index + k * rho->nx * rho->ny;
					  if (ChannelStopProfile == 0) // Square Well
						{
						  rho->data[index2] = ChannelStopCharge / CSChargeDepth;
						}
					  if (ChannelStopProfile == 1) // Gaussian
						{
						  rho->data[index2] = ChannelStopCharge / ChannelStopTotal / DeltaZ * (erf((rho->Z(rho->zp[k]+rho->dzp/2.0) - ChannelStopZmin) / (sqrt(2.0) * ChannelStopDepth)) - erf((rho->Z(rho->zp[k]-rho->dzp/2.0) - ChannelStopZmin) / (sqrt(2.0) * ChannelStopDepth)));
						}
					}
				}
			  else
				{
				  // This is the channel region
				  for (k=Channelkmin; k<Channelkmax+1; k++)
					{
					  DeltaZ = rho->Z(rho->zp[k] + rho->dzp / 2.0) - rho->Z(rho->zp[k] - rho->dzp / 2.0);
					  index2 = index + k * rho->nx * rho->ny;
					  if (ChannelProfile == 0) // Square Well
						{
						  rho->data[index2] = ChannelCharge / CChargeDepth;
						}
					  if (ChannelProfile == 1) // Gaussian
						{
						  rho->data[index2] = ChannelCharge / ChannelTotal / DeltaZ * (erf((rho->Z(rho->zp[k]+rho->dzp/2.0) - ChannelZmin) / (sqrt(2.0) * ChannelDepth)) - erf((rho->Z(rho->zp[k]-rho->dzp/2.0) - ChannelZmin) / (sqrt(2.0) * ChannelDepth)));
						}
					  /*
					  if (i == nx2 && j == ny2)
						{
						  //printf("DeltaZ = %.4f\n",DeltaZ);
						  TotalChannelCharge+=rho->data[index2] * DeltaZ;
						}Debug purposes */  
					}
				}
			}
		}
    }
  printf("Finished setting Fixed Charges, \n");
  //printf("TotalChannelCharge = %.6g\n", TotalChannelCharge);  
  fflush(stdout);
  return;
}

void MultiGrid::SetInitialElectrons(Array3D* rho, Array3D* elec)
{
  int i, j, k, n, q, index, index2;
  int CollectedChargeimin, CollectedChargeimax, CollectedChargeXWidth;
  CollectedChargeimin = rho->XIndex(CollectedChargeXmin);
  CollectedChargeimax = rho->XIndex(CollectedChargeXmax);
  CollectedChargeXWidth = CollectedChargeimax - CollectedChargeimin + 1;
  int CollectedChargejmin, CollectedChargejmax, CollectedChargeYWidth;
  CollectedChargejmin = rho->YIndex(CollectedChargeYmin);
  CollectedChargejmax = rho->YIndex(CollectedChargeYmax);
  CollectedChargeYWidth = CollectedChargejmax - CollectedChargejmin + 1;


  CollectedChargeimin = (int)(CollectedChargeXmin / GridSpacing);
  CollectedChargeimax = (int)(CollectedChargeXmax / GridSpacing);
  CollectedChargeXWidth = CollectedChargeimax - CollectedChargeimin + 1;
  CollectedChargejmin = (int)(CollectedChargeYmin / GridSpacing);
  CollectedChargejmax = (int)(CollectedChargeYmax / GridSpacing);
  CollectedChargeYWidth = CollectedChargejmax - CollectedChargejmin + 1;



  int CollectedChargekmin, CollectedChargekmax, CollectedChargeZWidth;
  CollectedChargekmin = rho->ZIndex(CollectedChargeZmin);
  CollectedChargekmax = rho->ZIndex(CollectedChargeZmax);
  CollectedChargeZWidth = CollectedChargekmax - CollectedChargekmin + 1;
  int PixX, PixY, FilledPixX, FilledPixY, Pixi, Pixj;
  double PixXmin, PixYmin, CollectCharge=0.0;
  printf("In SetMobileCharges, CollectedChargeimin = %d, CollectedChargeimax = %d CollectedChargeXWidth = %d grid cells\n",CollectedChargeimin, CollectedChargeimax, CollectedChargeXWidth);
  printf("In SetMobileCharges, CollectedChargejmin = %d, CollectedChargejmax = %d CollectedChargeYWidth = %d grid cells\n",CollectedChargejmin, CollectedChargejmax, CollectedChargeYWidth);
  printf("In SetMobileCharges, CollectedChargekmin = %d, CollectedChargekmax = %d CollectedChargeZWidth = %d grid cells\n",CollectedChargekmin, CollectedChargekmax, CollectedChargeZWidth);  
    // Potentials and Charges in Pixel Region
  double TotalElectrons = 0.0;
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<elec->nx; i++)
	{
	  if (elec->x[i] < PixelRegionLowerLeft[n][0] || elec->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<elec->ny; j++)
	    {
	      if (elec->y[j] < PixelRegionLowerLeft[n][1] || elec->y[j] > PixelRegionUpperRight[n][1])
		{ 
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * elec->nx;
	      PixX = (int)((elec->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)((elec->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
	      // First set the charges

	      for (q=0; q<NumberofFilledWells[n]; q++)
		{
		  CollectCharge = CollectedCharge[n][q] / ((double)CollectedChargeXWidth * (double)CollectedChargeYWidth * (double)CollectedChargeZWidth);		      
		  FilledPixX = (int)((FilledPixelCoords[n][q][0] - PixelRegionLowerLeft[n][0]) / PixelSize);
		  FilledPixY = (int)((FilledPixelCoords[n][q][1] - PixelRegionLowerLeft[n][1]) / PixelSize);
		  if (!(FilledPixX == PixX && FilledPixY == PixY))
		    {
		      continue;
		    } 
		  Pixi = (int)((elec->x[i] - PixXmin) / GridSpacing);
		  Pixj = (int)((elec->y[j] - PixYmin) / GridSpacing);
		  if (Pixi >= CollectedChargeimin && Pixi <= CollectedChargeimax && Pixj >= CollectedChargejmin && Pixj <= CollectedChargejmax)
		    {
		      for (k=CollectedChargekmin; k<CollectedChargekmax + 1; k++)
			{
			  index2 = index + k * elec->nx * elec->ny;
			  elec->data[index2] += CollectCharge;
			  TotalElectrons += CollectCharge;
			}
		    }
		}
	    }
	}
    }
  printf ("In SetInitialElectrons, CollectCharge = %.3f\n", CollectCharge);
  printf("Finished setting Initial Electrons, %.1f total electrons placed.\n", TotalElectrons);
  fflush(stdout);
  return;
}

void MultiGrid::SetInitialHoles(Array3D* rho, Array3D* hole)
{
  // This sets an initial guess at the free hole density by setting it to the
  // expected solution based on the 1D solution.
  // It is then adjusted by the routine AdjustHoles
  int i, j, k, n, index, index2;
  int PixX, PixY;
  double PixXmin, PixYmin, AddedHoles, TotalAddedHoles = 0.0;
  double ChargeDepth, ZL, Vdelta, Slope1, Slope2, SlopeDelta;
  double RhoChargeFactor = (QE*MICRON_PER_M / (EPSILON_0*EPSILON_SI)) / (GridSpacing * GridSpacing);
  ZL = rho->zmax - rho->Z(rho->zp[ChannelStopkmax] - rho->dzp / 2.0);
  Vdelta = BackgroundDoping / pow(MICRON_PER_CM, 3) * RhoChargeFactor * ZL * ZL * GridSpacing * GridSpacing;
  printf("In SetInitialHoles\n");

  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
	{
	  if (rho->x[i] < PixelRegionLowerLeft[n][0] || rho->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<rho->ny; j++)
	    {
	      if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
		{ 
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * rho->nx;
	      PixX = (int)((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)((rho->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
	      // Now set the charges
	      if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0)
		{
		  // This is the channel stop region
		  if (((rho->y[j] >= PixYmin + PixelSize/6.0 && rho->y[j] <= PixYmin + 5.0*PixelSize/6.0) && CollectingPhases == 2) || ((rho->y[j] >= PixYmin + PixelSize/3.0 && rho->y[j] <= PixYmin + 2.0*PixelSize/3.0) && CollectingPhases == 1))
		    {
		      // This is the Channel Stop Region under the collecting gates
		      Slope1 = (Vchannelstop - Vparallel_hi) / (rho->z[ChannelStopkmin] - rho->z[0]);
		    }
		  else
		    {
		      // This is the Channel Stop Region under the barrier gates
		      Slope1 = (Vchannelstop - Vparallel_lo) / (rho->z[ChannelStopkmin] - rho->z[0]);
		    }
		  for (k=ChannelStopkmin; k<ChannelStopkmax+1; k++)
		    {
		      index2 = index + k * rho->nx * rho->ny;
		      ChargeDepth = rho->Z(rho->zp[k] + rho->dzp / 2.0) - rho->Z(rho->zp[k] - rho->dzp / 2.0);  
		      SlopeDelta = rho->data[index2] * ChargeDepth * 2.0; // 2.0 Fudge factor not understood.
		      if (k == ChannelStopkmax)
			{
			  Slope1 = 0.0;
			  Slope2 = (Vbb - Vdelta) / ZL;
			  AddedHoles = (Slope1 - Slope2) / RhoChargeFactor;
			  AddedHoles += -rho->data[index2] * ChargeDepth / RhoChargeFactor; // Add in charge to compensate fixed charge
			  AddedHoles = max(AddedHoles, -hole->data[index2]);			  
			  hole->data[index2] = AddedHoles;
			  TotalAddedHoles += AddedHoles;
			}
		      else if ((Slope1 - SlopeDelta) < 0.0)
			{
			  //  Not enough fixed charge yet
			  Slope1 = Slope1 - SlopeDelta;
			  hole->data[index2] = 0.0;
			}
		      else
			{
			  // Now have enough fixed charge
			  Slope2 = 0.0;
			  AddedHoles = (Slope1 - Slope2) / RhoChargeFactor;
			  AddedHoles += -rho->data[index2] * ChargeDepth / RhoChargeFactor; // Add in charge to compensate fixed charge
			  AddedHoles = max(AddedHoles, -hole->data[index2]);			  
			  hole->data[index2] = AddedHoles;
			  TotalAddedHoles += AddedHoles;
			  Slope1 = 0.0;
			}
		    } // ends k
		} // ends channel stop region
	      else
		{
		  // This is the channel region
		  continue;
		}
	    } // ends j
	} // ends i
    } // ends n
  double TotHoles = 0.0;
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
	{
	  for (k=0; k<hole->nz; k++)
	    {
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      TotHoles += hole->data[index];
	    }
	}
    }
  printf("Finished setting Initial holes, %.6g total holes placed.\n", TotalAddedHoles);
  fflush(stdout);
  return;

}

void MultiGrid::AdjustHoles(Array3D* phi, Array3D* rho, Array3D* hole)
{
  // This routine adjusts the free hole density to "pin" the potential in
  // the channel stop region to be equal to Vchannelstop.
  // It needs to be run iteratively to converge to a solution.
  int i, j, k, n, index, index2;
  int PixX, PixY;
  double PixXmin, PixYmin, AddedHoles, MinusHoles, PartialHoles, TotalAddedHoles = 0.0, ChargeIncrement;
  //ChargeIncrement = -5.0 * ChannelStopCharge / CSChargeDepth / (double)(pow(ScaleFactor, 3.0));
  ChargeIncrement = 20.0 / (double)(pow(ScaleFactor, 3.0));  
  printf("ChargeIncrement = %f\n",ChargeIncrement);
  //  This sets how rapidly we converge to a solution.
  double MaxIncrement = 10000.0;
  double TotHoles = 0.0;
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
	{
	  for (k=0; k<hole->nz; k++)
	    {
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      TotHoles += hole->data[index];
	    }
	}
    }
  printf("Starting adjusting Mobile holes, %.6g total holes.\n", TotHoles);
  
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
	{
	  if (rho->x[i] < PixelRegionLowerLeft[n][0] || rho->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<rho->ny; j++)
	    {
	      if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
		{ 
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * rho->nx;
	      PixX = (int)((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)((rho->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
	      // Now set the charges
	      if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0)
		{
		  // This is the Channel Stop Region
		  for (k=ChannelStopkmin; k < ChannelStopkmax+1; k++)
		    {
		      index2 = index + k * rho->nx * rho->ny;		      
		      if (phi->data[index2] < Vchannelstop)			
			{
			  // Region needs more free holes
			  if ((Vchannelstop - phi->data[index2]) > 0.5)
			    {
			      AddedHoles = ChargeIncrement * (Vchannelstop - phi->data[index2]);
			    }
			  else
			    {
			      AddedHoles = 2.0 * ChargeIncrement * (Vchannelstop - phi->data[index2]);
			    }
			    AddedHoles = min(MaxIncrement, AddedHoles);
			  hole->data[index2] += AddedHoles;			  
			  TotalAddedHoles += AddedHoles;
			}
		      else
			{
			  //phi->data[index2] > Vchannelstop
			  // Region has too many free holes
			  if ((phi->data[index2] - Vchannelstop) > 0.5)
			    {
			      MinusHoles = ChargeIncrement * (phi->data[index2] - Vchannelstop);			  			    }
			  else
			    {
			      MinusHoles = 2.0 * ChargeIncrement * (phi->data[index2] - Vchannelstop);			  			    }
			  MinusHoles = min(MaxIncrement, MinusHoles);
			  if (MinusHoles < hole->data[index2])
			    {
			      // We're less than the hole value, so just subtract it.
			      hole->data[index2] -= MinusHoles;			  
			      TotalAddedHoles -= MinusHoles;
			    }
			  else
			    {
			      //  We're over the value to bring this to zero, so bring this cell to zero
			      PartialHoles = hole->data[index2];
			      hole->data[index2] -= PartialHoles;
			      TotalAddedHoles -= PartialHoles;			      
			    }
			}// ends else
		    }// ends k
		}// ends channel stop region
	      else
		{
		  // This is the channel region
		  continue;
		}
	    }
	}
    }
  TotHoles = 0.0;
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
	{
	  for (k=0; k<hole->nz; k++)
	    {
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      TotHoles += hole->data[index];
	    }
	}
    }
  printf("Finished adjusting Mobile Holes, %.6g added holes, %.6g total holes.\n", TotalAddedHoles, TotHoles);
  fflush(stdout);
  return;
}

void MultiGrid::SOR(Array3D* phi, Array3D* rho, Array2D* BCType, double w)
{
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  double newphi;
  double omw, w6, hsquared;
  w6 = w / 6.0;
  int i, j, k, im, ip, j0, jm, jp, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  omw = 1.0 - w;
  for (i=0; i<phi->nx; i++)
    {
      if (XBCType == 0)
	{
	  if (i == 0) {im = 0;} else {im = 1;} // Free BC
	  if (i == phi->nx-1) {ip = 0;} else {ip = 1;} // Free BC
	}
      else
	{
	  if (i == 0) {im = -phi->nx + 2;} else {im = 1;} // Periodic BC
	  if (i == phi->nx-1) {ip = -phi->nx + 2;} else {ip = 1;} // Periodic BC
	}
      for (j=0; j<phi->ny; j++)
	{
	  j0 = j * phi->nx;
      if (YBCType == 0)
	{
	  if (j == 0) {jm = 0;} else {jm = phi->nx;} // Free BC
	  if (j == phi->ny-1) {jp = 0;} else {jp = phi->nx;} // Free BC
	}
      else
	{
	  if (j == 0) {jm = (-phi->ny + 2) * phi->nx;} else {jm = phi->nx;} // Periodic BC
	  if (j == phi->ny-1) {jp = (-phi->ny + 2) * phi->nx;} else {jp = phi->nx;} // Periodic BC
	}
	  ind = i + j0;
	  if (BCType->data[ind] == 1) // Free BC at z = 0
	    {
	      indmx = ind - im;
	      indpx = ind + ip;
	      indmy = ind - jm;
	      indpy = ind + jp;
	      indmz = ind;
	      indpz = ind + nxy;
	      newphi = omw * phi->data[ind] + w6 * (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->data[indmz]+phi->data[indpz] + hsquared * rho->data[ind]);
	      phi->data[ind] = newphi;
	    }
	  for (k=1; k<phi->nz-1; k++)
	    {
	      w6 = w / (4.0 + phi->zplus[k] + phi->zminus[k]);
	      ind = i + j0 + k * nxy;
	      indmx = ind - im;
	      indpx = ind + ip;
	      indmy = ind - jm;
	      indpy = ind + jp;
	      indmz = ind - nxy;
	      indpz = ind + nxy;
	      newphi = omw * phi->data[ind] + w6 * (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->zminus[k]*phi->data[indmz]+phi->zplus[k]*phi->data[indpz] + hsquared * rho->data[ind]);
	      phi->data[ind] = newphi;
	    }
	}
    }
  return;
}


double MultiGrid::Error(Array3D* phi, Array3D* rho)
{
  double newphi, error = 0.0;
  double hsquared;
  int i, j, k, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  for (i=1; i<phi->nx-1; i++)
    { 
      for (j=1; j<phi->ny-1; j++)
	{
	  for (k=1; k<phi->nz-1; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      indmx = ind - 1;
	      indpx = ind + 1;
	      indmy = ind - phi->nx;
	      indpy = ind + phi->nx;
	      indmz = ind - nxy;
	      indpz = ind + nxy;
	      newphi = (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->zminus[k]*phi->data[indmz]+phi->zplus[k]*phi->data[indpz] + hsquared * rho->data[ind]) / (4.0+phi->zplus[k]+phi->zminus[k]);
	      error = max(error, fabs(phi->data[ind] - newphi));
	    }
	}
    }
  //printf("Grid Size = %d, error = %.3g\n",phi->nx, error);
  return error;
}


void MultiGrid::Restrict(Array3D* phi, Array3D* newphi, Array3D* rho, Array3D* newrho, Array2D* BCType, Array2D* newBCType)
{
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  int i, j, k, im, i0, ip, jm, j0, jp, km, k0, kp, nxy, newnxy;
  int newindex, ind000;
  int ind00p, ind0p0, indp00;
  int ind00m, ind0m0, indm00;
  int ind0mm, ind0mp, ind0pm, ind0pp;
  int indm0m, indm0p, indp0m, indp0p;
  int indmm0, indmp0, indpm0, indpp0;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  double phisum, rhosum;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;

  for (i=0; i<newphi->nx; i++)
    {
      if (XBCType == 0) // Free BC
	{
	  if (i == 0) {i0 = 0; im = 0; ip = 1;}
	  else if (i == newphi->nx - 1) {i0 = phi->nx - 1; im = phi->nx - 2; ip = phi->nx - 1;}
	  else {i0 = 2 * i; im = 2 * i - 1; ip = 2 * i + 1;}
	}
      else // Periodic BC
	{
	  if (i == 0) {i0 = 0; im = phi->nx - 2; ip = 1;}
	  else if (i == newphi->nx - 1) {i0 = phi->nx - 1; im = phi->nx - 2; ip = 1;}
	  else {i0 = 2 * i; im = 2 * i - 1; ip = 2 * i + 1;}
	}
      for (j=0; j<newphi->ny; j++)
	{
	  if (YBCType == 0) // Free BC
	    {
	      if (j == 0) {j0 = 0; jm = 0; jp = 1;}
	      else if (j == newphi->ny - 1) {j0 = phi->ny - 1; jm = phi->ny - 2; jp = phi->ny - 1;}
	      else {j0 = 2 * j; jm = 2 * j - 1; jp = 2 * j + 1;}
	    }
	  else // Periodic BC
	    {
	      if (j == 0) {j0 = 0; jm = phi->ny - 2; jp = 1;}
	      else if (j == newphi->ny - 1) {j0 = phi->ny - 1; jm = phi->ny - 2; jp = 1;}
	      else {j0 = 2 * j; jm = 2 * j - 1; jp = 2 * j + 1;}
	    }
	  newBCType->data[i + j * newphi->nx] = BCType->data[i0 + j0 * phi->nx];
	  for (k=0; k<newphi->nz; k++)
	    {
	      if (k == 0) {k0 = 0; km = 0; kp = 0;}
	      else if (k == newphi->nz - 1) {k0 = phi->nz - 1; km = phi->nz - 1; kp = phi->nz - 1;}
	      else {k0 = 2 * k; km = 2 * k - 1; kp = 2 * k + 1;}
	      newindex = i + j*newphi->nx + k*newnxy;
	      ind000 = i0 + j0 * phi->nx + k0 * nxy;
	      phisum = 8.0 * phi->data[ind000];
	      rhosum = 8.0 * rho->data[ind000];
	      indm00 = im + j0 * phi->nx + k0 * nxy;
	      indp00 = ip + j0 * phi->nx + k0 * nxy;
	      ind0m0 = i0 + jm * phi->nx + k0 * nxy;
	      ind0p0 = i0 + jp * phi->nx + k0 * nxy;
	      ind00m = i0 + j0 * phi->nx + km * nxy;
	      ind00p = i0 + j0 * phi->nx + kp * nxy;
	      phisum += 4.0 * (phi->data[indm00] + phi->data[indp00] + phi->data[ind0m0] + phi->data[ind0p0] + phi->data[ind00m] + phi->data[ind00p]);
	      rhosum += 4.0 * (rho->data[indm00] + rho->data[indp00] + rho->data[ind0m0] + rho->data[ind0p0] + rho->data[ind00m] + rho->data[ind00p]);
	      indmm0 = im + jm * phi->nx + k0 * nxy;
	      indmp0 = im + jp * phi->nx + k0 * nxy;
	      indpm0 = ip + jm * phi->nx + k0 * nxy;
	      indpp0 = ip + jp * phi->nx + k0 * nxy;
	      phisum += 2.0 * (phi->data[indmm0] + phi->data[indmp0] + phi->data[indpm0] + phi->data[indpp0]);
	      rhosum += 2.0 * (rho->data[indmm0] + rho->data[indmp0] + rho->data[indpm0] + rho->data[indpp0]);
	      indm0m = im + j0 * phi->nx + km * nxy;
	      indm0p = im + j0 * phi->nx + kp * nxy;
	      indp0m = ip + j0 * phi->nx + km * nxy;
	      indp0p = ip + j0 * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[indm0m] + phi->data[indm0p] + phi->data[indp0m] + phi->data[indp0p]);
	      rhosum += 2.0 * (rho->data[indm0m] + rho->data[indm0p] + rho->data[indp0m] + rho->data[indp0p]);
	      ind0mm = i0 + jm * phi->nx + km * nxy;
	      ind0mp = i0 + jm * phi->nx + kp * nxy;
	      ind0pm = i0 + jp * phi->nx + km * nxy;
	      ind0pp = i0 + jp * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[ind0mm] + phi->data[ind0mp] + phi->data[ind0pm] + phi->data[ind0pp]);
	      rhosum += 2.0 * (rho->data[ind0mm] + rho->data[ind0mp] + rho->data[ind0pm] + rho->data[ind0pp]);
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      phisum += phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp];
	      rhosum += rho->data[indmmm] + rho->data[indmmp] + rho->data[indmpm] + rho->data[indmpp] + rho->data[indpmm] + rho->data[indpmp] + rho->data[indppm] + rho->data[indppp];
	      newphi->data[newindex] = phisum / 64.0;
	      newrho->data[newindex] = rhosum / 64.0;
	    }
	}
    }
  return;
}


void MultiGrid::Prolongate(Array3D* phi, Array3D* newphi, Array2D* newBCType)
{
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  int i, j, k, im, ip, jm, jp, km, kp, nxy, newnxy;
  int newindex;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;
  for (i=0; i<newphi->nx; i++)
    {
      im = max(0,rint((float)i / 2.0 - 0.1));
      ip = min(newphi->nx-1,rint((float)i / 2.0 + 0.1));
      for (j=0; j<newphi->ny; j++)
	{
	  jm = max(0,rint((float)j / 2.0 - 0.1));
	  jp = min(newphi->ny-1,rint((float)j / 2.0 + 0.1));
	  if (newBCType->data[i + j * newphi->nx] == 1) // Free BC at z = 0
	    {
	      newindex = i + j * newphi->nx;
	      indmmm = im + jm * phi->nx;
	      indmmp = im + jm * phi->nx;
	      indmpm = im + jp * phi->nx;
	      indmpp = im + jp * phi->nx;
	      indpmm = ip + jm * phi->nx;
	      indpmp = ip + jm * phi->nx;
	      indppm = ip + jp * phi->nx;
	      indppp = ip + jp * phi->nx;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	  for (k=1; k<newphi->nz-1; k++)
	    {
	      km = rint((float)k / 2.0 - 0.1);
	      kp = rint((float)k / 2.0 + 0.1);
	      newindex = i + j * newphi->nx + k * newnxy;
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	}
    }
  return;
}

void MultiGrid::VCycle(Array3D** phi, Array3D** rho, Array2D** BCType, double w, int nsteps, int ncycle)
{
  // Iterates a few steps (4) at each grid on the way down to smooth, then
  // iterates the coarsest grid down to machine precision,
  // then does a given number of steps (ncycle) at each finer scale on the way up.
  int i, j, niter;
  double error;
  
  for (i=0; i<nsteps; i++)
    {
      Restrict(phi[i], phi[i+1], rho[i], rho[i+1], BCType[i], BCType[i+1]);
      for (j=0; j<4; j++)
	{
	  SOR(phi[i+1], rho[i+1], BCType[i+1], w);
	}
    }
  error = 100.0; niter = 0;
  while (error > 1.0E-12)
    {
      niter += 1;
      SOR(phi[nsteps], rho[nsteps], BCType[nsteps], w);
      error = Error(phi[nsteps], rho[nsteps]);
    }
  printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g\n",phi[nsteps]->nx-1,phi[nsteps]->ny-1,phi[nsteps]->nz-1,niter,error);
  fflush(stdout);
  for (i=nsteps; i>0; i--)
    {
      Prolongate(phi[i], phi[i-1], BCType[i-1]);
      niter = ncycle * (int)pow(2,i-1);
      for (j=0; j<niter; j++)
	{
	  SOR(phi[i-1], rho[i-1], BCType[i-1], w);
	}
      error = Error(phi[i-1], rho[i-1]);
      printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g\n",phi[i-1]->nx-1,phi[i-1]->ny-1,phi[i-1]->nz-1,niter,error);
      fflush(stdout);
    }
  return;
}


void MultiGrid::WriteOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This writes the data to the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  int* int_attr_data  = new int[3];
  double* double_attr_data  = new double[3];
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;
  
  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nz; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      flipped_data[flipped_index] = array->data[index];
	    }
	}
    }
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname;
  WriteHDF5File3(filename, hdfname, array->nx, array->ny, array->nz, flipped_data);
  // Now we write the attributes
  int_attr_data[0] = array->nx; int_attr_data[1] = array->ny; int_attr_data[2] = array->nz;
  WriteHDF5IntAttribute(filename, hdfname, "Dimension", 3, int_attr_data);
  double_attr_data[0] = array->xmin; double_attr_data[1] = array->ymin; double_attr_data[2] = array->zmin;
  WriteHDF5DoubleAttribute(filename, hdfname, "Lower_Left", 3, double_attr_data);
  double_attr_data[0] = array->xmax; double_attr_data[1] = array->ymax; double_attr_data[2] = array->zmax;
  WriteHDF5DoubleAttribute(filename, hdfname, "Upper_Right", 3, double_attr_data);

  delete[] flipped_data;
  return;
}

void MultiGrid::Gradient(Array3D* phi, Array3D** E)
{
  int i, j, k, nxy, ind;
  nxy = phi->nx * phi->ny;

  // Ex 
  for (j=0; j<phi->ny; j++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind]) / phi->dx;
	  ind = 1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-2 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind] - phi->data[ind-1]) / phi->dx;
	  for (i=2; i<phi->nx-2; i++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[0]->data[ind] = (-phi->data[ind+2] + 8.0 * phi->data[ind+1] - 8.0 * phi->data[ind-1] + phi->data[ind-2]) / (12.0 * phi->dx);
	    }
	}
    }

    // Ey 
  for (i=0; i<phi->nx; i++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = i + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind]) / phi->dy;
	  ind = i + phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-2) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-1) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind] - phi->data[ind-phi->nx]) / phi->dy;
	  for (j=2; j<phi->ny-2; j++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[1]->data[ind] = (-phi->data[ind+2*phi->nx] + 8.0 * phi->data[ind+phi->nx] - 8.0 * phi->data[ind-phi->nx] + phi->data[ind-2*phi->nx]) / (12.0 * phi->dy);
	    }
	}
    }

  // Ez 
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  ind = i + j * phi->nx;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind]) / phi->dzp * phi->dzpdz[0];
	  ind = i + j * phi->nx + nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp)  * phi->dzpdz[1];
	  ind = i + j * phi->nx + (phi->nz-2) * nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp) * phi->dzpdz[phi->nz-2];
	  ind = i + j * phi->nx + (phi->nz-1) * nxy;
	  E[2]->data[ind] = (phi->data[ind] - phi->data[ind-nxy]) / phi->dzp * phi->dzpdz[phi->nz-1];
	  for (k=2; k<phi->nz-2; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[2]->data[ind] = (-phi->data[ind+2*nxy] + 8.0 * phi->data[ind+nxy] - 8.0 * phi->data[ind-nxy] + phi->data[ind-2*nxy]) / (12.0 * phi->dzp) * phi->dzpdz[k];
	    }
	}
    }
  return;
}

void MultiGrid::Trace(double* point, int bottomsteps, bool savecharge, double bottomcharge, ofstream& file)
{
  // This traces an electron down to the bottom, saving path info if requested
  // Diffusion has now been added. This version recalculates mu at each point.
  // And iterates bottomsteps steps after reaching the bottom.
  // If savecharge is true, it finds and stores the self-consistent charge locations 
  int i, j, k, nsteps = 0, nstepsmax = 10000;
  bool ReachedBottom = false;
  double mu, E2, Emag, ve, vth, tau, Tscatt;
  double theta, phiangle, zmin, zbottom;
  double x, y, z;
  zmin = (E[0]->z[Channelkmax] + E[0]->z[Channelkmin]) / 2.0;
  zbottom = E[0]->Z(E[0]->zp[Channelkmin] - E[0]->dzp / 2.0 + 0.01);
  x = point[0]; y = point[1]; z = point[2];
  double*  E_interp = new double[3];
  
  while (nsteps < nstepsmax)    
    {
      nsteps += 1;
      E2 = 0.0;
      for (i=0; i<3; i++)
	{
	  E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
	  E2 += E_interp[i] * E_interp[i];
	}
      Emag = max(0.1, sqrt(E2));
      mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
      ve = mu * MICRON_PER_CM * MICRON_PER_CM; // Drift Velocity Factor (Velocity / E)
      vth = sqrt(3.0 * KBOLTZMANN * CCDTemperature / ME)  * MICRON_PER_M * DiffMultiplier; // Thermal Velocity
      vth = vth / sqrt((double)NumDiffSteps);
      tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time

      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
      point[0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
      point[1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
      point[2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;      
      if (LogPixelPaths == 1)
	{
	  file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
	}
      if (point[2] < zmin && !ReachedBottom)
	{
	  ReachedBottom = true;
	  nstepsmax = nsteps + bottomsteps + bottomsteps / 10; 
	  // After reaching bottom, iterate bottomsteps (*1.1) more steps.
	  // The first bottomsteps/10 steps are to let it settle to an
	  // equilibrium location, then we start logging the charge location
	}
      if (ReachedBottom && nsteps > nstepsmax - bottomsteps)
	{
	  //  Start logging location after bottomsteps / 10.
	  point[2] = max(zbottom, point[2]);
	  i = E[0]->XIndex(point[0]);
	  j = E[0]->YIndex(point[1]);
	  k = E[0]->ZIndex(point[2]);
	  if (i > 0 && i < elec[0]->nx-1 && j > 0 && j < elec[0]->ny-1 && k < elec[0]->nz-1 && savecharge)
	    {
	      elec[0]->data[i + j * elec[0]->nx + k * elec[0]->nx * elec[0]->ny] += bottomcharge;// Add bottomcharge to this grid cell
	    }
	}
    }
  delete[] E_interp;
  return;
}

void MultiGrid::TraceSpot(int m)
{
  // This builds up a Gaussian spot with given center (Xoffset, Yoffset) and SigmaX and SigmaY
  double x, y, z, rsq, v1, v2, fac, xwindow, ywindow, xcenter, ycenter;
  int n;
  int bottomsteps = 1000;
  double bottomcharge = .001;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);      
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name);
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  file  << setw(15) << "xin" << setw(15) << "yin" << setw(15) << "zin" << setw(15) << "xout" << setw(15) << "yout" << setw(15) << "zout" << endl;
  xwindow = PixelBoundaryUpperRight[0] - PixelBoundaryLowerLeft[0];
  ywindow = PixelBoundaryUpperRight[1] - PixelBoundaryLowerLeft[1];
  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;  
  for (n=0; n<NumElec; n++)
    {
      //  Use Box-Muller algorithm to generate two Gaussian random numbers
      rsq = 1000.0;      
      while (rsq >= 1.0 || rsq == 0.0)
	{
	  v1 = 2.0 * drand48() - 1.0;
	  v2 = 2.0 * drand48() - 1.0;
	  rsq = v1*v1 + v2 *v2;
	}
      fac = sqrt(-2.0 * log(rsq) / rsq);
      x = xcenter + Sigmax * v1 * fac;
      y = ycenter + Sigmay * v2 * fac;
      point[0] = x;
      point[1] = y;
      z = ElectronZ0Fill;
      point[2] = z;
      Trace(point, bottomsteps, true, bottomcharge, file);
      // Trace returns the final location
      file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  return;
}

void MultiGrid::TraceMultipleSpots(int m)
{
  // This builds up multiple Gaussian spots
  // Extents still hard coded - needs some work
  double x, y, z, rsq, v1, v2, fac, xcenter, ycenter;
  int i, j, n, Nume;
  int bottomsteps = 1000;
  double bottomcharge = .001;
  double* point = new double[3];
  ofstream dummyfile;
  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;  
  for (i=-4; i<5; i++)
    {
      for (j=-4; j<5; j++)
	{
	  if (i == 0 && j == 0)
	    {
	      Nume = NumElec + NumElec / 10;
	    }
	  else
	    {
	      Nume = NumElec;
	    }
	  for (n=0; n<Nume; n++)
	    {
	      //  Use Box-Muller algorithm to generate two Gaussian random numbers
	      rsq = 1000.0;      
	      while (rsq >= 1.0 || rsq == 0.0)
		{
		  v1 = 2.0 * drand48() - 1.0;
		  v2 = 2.0 * drand48() - 1.0;
		  rsq = v1*v1 + v2 *v2;
		}
	      fac = sqrt(-2.0 * log(rsq) / rsq);
	      x = xcenter + (double)i * PixelSize + Sigmax * v1 * fac;
	      y = ycenter + (double)j * PixelSize + Sigmay * v2 * fac;
	      point[0] = x;
	      point[1] = y;
	      z = ElectronZ0Fill;
	      point[2] = z;
	      Trace(point, bottomsteps, true, bottomcharge, dummyfile);
	      // Trace returns the final location
	    }
	}
    }
  delete[] point;
  return;
}

void MultiGrid::TraceGrid(int m)
{
  // This traces a grid of starting electron locations.
  double x, y, z;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);      
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name);
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  file  << setw(15) << "xin" << setw(15) << "yin" << setw(15) << "zin" << setw(15) << "xout" << setw(15) << "yout" << setw(15) << "zout" << endl;

  x = PixelBoundaryLowerLeft[0] + PixelBoundaryStepSize[0] / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      y = PixelBoundaryLowerLeft[1] + PixelBoundaryStepSize[1] / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  point[0] = x;
	  point[1] = y;
	  z = ElectronZ0Fill;
	  point[2] = z;
	  Trace(point, 100, false, 0.0, file);
	  file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
	  y += PixelBoundaryStepSize[1];
	}
      x += PixelBoundaryStepSize[0];
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  return;
}

void MultiGrid::TraceRegion(int m)
{
  // This traces a random set of starting electron locations within the PixelBoundary.
  double x, y, z, boxx, boxy;
  int n;
  boxx = PixelBoundaryUpperRight[0] - PixelBoundaryLowerLeft[0];
  boxy = PixelBoundaryUpperRight[1] - PixelBoundaryLowerLeft[1];  
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);      
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name);
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  file  << setw(15) << "xin" << setw(15) << "yin" << setw(15) << "zin" << setw(15) << "xout" << setw(15) << "yout" << setw(15) << "zout" << endl;
  for (n=0; n<NumElec; n++)
    {
      if (n%1000==0)
	{
	  printf("Finished %d electrons\n",n);
	  fflush(stdout);
	}
      x = PixelBoundaryLowerLeft[0] + drand48() * boxx;
      y = PixelBoundaryLowerLeft[1] + drand48() * boxy;
      point[0] = x;
      point[1] = y;
      z = ElectronZ0Fill;
      point[2] = z;
      Trace(point, 100, false, 0.0, file);
      file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  return;
}

void MultiGrid::FindEdge(double* point, double theta, ofstream& file)
{
  // This finds the edge of the pixel through binary search given a starting point and a line angle
  int nsteps, pixx, pixy, lastpixx, lastpixy, newpixx, newpixy;
  double sinth, costh, x, y, z0, deltar, tolerance, sign;
  sinth = sin(theta);
  costh = cos(theta);
  z0 = point[2];
  x = point[0]; y = point[1];
  pixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSize);
  pixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSize);
  lastpixx = pixx; lastpixy = pixy;
  deltar = 1.0;
  tolerance = 0.0001;
  nsteps = 0;
  while (fabs(deltar) > tolerance)
    {
      nsteps += 1;
      if (nsteps >200)
	{
	  printf("Too many steps in edge finding, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, nsteps,x,y);
	  fflush(stdout);
	  break;
	}
      x += deltar * costh;
      y += deltar * sinth;
      point[0] = x;
      point[1] = y;
      point[2] = z0;
      Trace(point, 10, false, 0.0, file);      
      newpixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSize);
      newpixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSize);
      //printf("Finding edge, newpixx = %d, newpixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",newpixx, newpixy, theta, nsteps,point[0],point[1]);
      if (newpixx == pixx && newpixy == pixy)
	{
	  sign = 1.0;
	}
      else
	{
	  sign = -1.0;
	}
      if (newpixx == lastpixx && newpixy == lastpixy)
	{
	  deltar = sign * fabs(deltar);
	}
      else
	{
	  deltar = sign * fabs(deltar) / 2.0;
	}
      lastpixx = newpixx; lastpixy = newpixy;
    }
  //printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, nsteps,x,y);
  //fflush(stdout);
  point[0] = x;
  point[1] = y;
  point[2] = z0;
  return;
}

void MultiGrid::FindCorner(double* point, double* theta, ofstream& file)
{
  // This finds the corner of the pixel through binary search given a starting point and a line angle
  int nsteps;
  double theta0, delta_theta, r, rm, r0, rp, deltar, tolerance, x0, y0, z0;
  theta0 = *theta;
  delta_theta = 0.10;
  z0 = point[2];
  x0 = point[0]; y0 = point[1];
  deltar = 1.0;
  tolerance = 0.0001;
  nsteps = 0;
  while (fabs(deltar) > tolerance)
    {
      nsteps += 1;
      if (nsteps >200)
	{
	  printf("Too many steps in corner finding, theta = %.3f, %d steps, x0 = %.3f, y0 = %.3f\n",theta0, nsteps,x0,y0);
	  fflush(stdout);
	  break;
	}
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      r0 = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 + delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rp = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 - delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rm = r;
      //printf("Step = %d, theta = %.5f, r0 = %.4f, rm = %.4f, rp = %.4f\n",nsteps,theta,r0,rm,rp);
      if (rp > r0)
	{
	  theta0 = theta0 + delta_theta;
	  deltar = fabs(r0 - rp);
	}
      else if (rm > r0)
	{
	  theta0 = theta0 - delta_theta;
	  deltar = fabs(r0 - rm);
	}
      else if (rp > rm)
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rp);
	}
      else
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rm);
	}
    }
  //printf("Found corner, theta = %.4f, x = %.4f, y = %.4f, after %d steps\n",theta0, point[0], point[1], nsteps);
  //fflush(stdout);
  *theta = theta0;
  return ;
}

void MultiGrid::CalculatePixelAreas(int m)
{
  // This finds the pixel vertices and the areas of the pixel grid.
  int OldLogPixelPaths, k, n, pixx, pixy;
  //Turn off Pixel Path logging
  OldLogPixelPaths = LogPixelPaths;
  LogPixelPaths = 0;
  //Turn off Diffusion
  double OldDiffMultiplier;
  OldDiffMultiplier = DiffMultiplier;
  DiffMultiplier = 0.0;
  double x, y, xb, yb, theta, theta0, theta1, dtheta, area;
  double* point = new double[3];
  string dummy = "dummy";
  string ptsfilename = (dummy);
  ofstream ptsfile; //Not needed, but we need to pass something to the Trace subroutine

  // Now calculate the pixel vertices
  Polygon** polyarray = new Polygon*[PixelBoundaryNx * PixelBoundaryNy];
  x = PixelBoundaryLowerLeft[0] + PixelSize / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSize);
      y = PixelBoundaryLowerLeft[1] + PixelSize / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSize);
	  polyarray[pixx + PixelBoundaryNx * pixy] = new Polygon(4 * NumVertices + 4);
	  //First, find the four corners
	  for (n=1; n<8; n+=2)
	    {
	      theta = (double)n * pi / 4.0;
	      point[0] = x;
	      point[1] = y;
	      point[2] = ElectronZ0Area;
	      FindCorner(point, &theta, ptsfile);
	      Point* two_d_point = new Point(point[0], point[1], theta);	      
	      polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
	    }
	  // Now find NumVertices points along each edge
	  for (n=0; n<4; n++) // Four corners
	    {
	      for (k=0; k<NumVertices; k++)
		{
		  if (n == 3) // Need this to make angles continuous through zero
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta - 2.0 * pi;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->theta;		  
		    }
		  else
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[(n+1)]->theta;		  
		    }
		  dtheta = (theta1 - theta0) / ((double)NumVertices + 1.0);
		  theta  = theta0  + ((double)k + 1.0) * dtheta;
		  point[0] = x;
		  point[1] = y;
		  point[2] = ElectronZ0Area;
		  FindEdge(point, theta, ptsfile);
		  //printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, x = %.3f, y = %.3f\n",pixx, pixy, theta, point[0],point[1]);
		  Point* two_d_point = new Point(point[0], point[1], theta);	      
		  polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
		}
	    }
	  printf("Finished vertex finding for pixel, pixx = %d, pixy = %d\n",pixx, pixy);
	  printf("Corners at (%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f)\n",
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->y, 
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->y, 
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->y, 
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->y); 
	  fflush(stdout);
	  y += PixelSize;
	}
      x += PixelSize;
    }

  // Now calculate and print out the pixel areas
  string underscore = "_", slash = "/", vertexname = "Vertices", areaname = "Area";
  string StepNum = boost::lexical_cast<std::string>(m);      
  string areafilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+areaname);  
  ofstream areafile;  
  areafile.open(areafilename.c_str());
  areafile.setf(ios::fixed);
  areafile.setf(ios::showpoint);
  areafile.setf(ios::left);
  areafile.precision(4);
  areafile  << setw(15) << "Nx" << setw(15) << "Ny" << setw(15) << "Area" << endl;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  area = polyarray[pixx + PixelBoundaryNx * pixy]->Area();	  
	  printf("Found area, pixx = %d, pixy = %d, area = %.3f\n",pixx, pixy, area);
	  fflush(stdout);
	  areafile  << setw(15) << pixx << setw(15) << pixy << setw(15) << area << endl;
	}
    }
  areafile.close();
  printf("Finished writing grid file - %s\n",areafilename.c_str());

  // Now print out the pixel vertices
  string vertexfilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+vertexname);  
  ofstream vertexfile;  
  vertexfile.open(vertexfilename.c_str());
  vertexfile.setf(ios::fixed);
  vertexfile.setf(ios::showpoint);
  vertexfile.setf(ios::left);
  vertexfile.precision(4);
  vertexfile  << setw(15) << "X0" << setw(15) << "Y0" << setw(15)<< "Theta" << setw(15) << "X" << setw(15) << "Y" << endl;
  x = PixelBoundaryLowerLeft[0] + PixelSize / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSize);
      y = PixelBoundaryLowerLeft[1] + PixelSize / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSize);
	  for (n=0; n<polyarray[pixx + PixelBoundaryNx * pixy]->npoints; n++)
	    {
	      xb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->x;
	      yb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->y;
	      theta = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;	      	      
	      vertexfile  << setw(15) << x << setw(15) << y << setw(15)<< theta << setw(15)<< xb << setw(15) << yb << endl;
	    }
	  y += PixelSize;
	}
      x += PixelSize;
    }
  vertexfile.close();  
  ptsfile.close();
  printf("Finished writing grid file - %s\n",vertexfilename.c_str());
  fflush(stdout);
  // Clean up
  delete[] point;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  delete polyarray[pixx + PixelBoundaryNx * pixy];
	}
    }
  delete[] polyarray;
  //Put Pixel path logging back where it was
  LogPixelPaths = OldLogPixelPaths;
  //Turn Diffusion back on
  DiffMultiplier = OldDiffMultiplier;
  return;
}

void MultiGrid::AddDipolePotentials(Array3D* phi)
{
  int i, j, k, ndipole, nimage, index;
  double r2, r3, deltaphi, p, px, py, pz, dipolesign;
  double DipoleFactor =  (QE/(4.0 * pi * EPSILON_0*EPSILON_SI)) * MICRON_PER_M;
  // DipoleFactor converts (charge in e-) * (separation in microns) into the appropriate units
  for (i=0; i<phi->nx; i++)
    { 
      for (j=0; j<phi->ny; j++)
	{
	  for (k=1; k<phi->nz-1; k++)
	    {
	      index = i + j * phi->nx + k * phi->nx * phi->ny;
	      for (ndipole=0; ndipole<NumberofDipoles; ndipole++)
		{
		  p = DipoleFactor * (double)DipoleCharge[ndipole] * (2.0 * DipoleZLocation[ndipole]);
		  // p is the dipole strength
		  // p(x,y,z) are the vector coordinates from the dipole to the grid point 
		  px = phi->x[i] - DipoleCoords[ndipole][0];
		  py = phi->y[j] - DipoleCoords[ndipole][1];
		  r2 = px*px + py*py;
		  for (nimage= - NumberofDipoleImages; nimage<NumberofDipoleImages+1; nimage++)
		    {
		      if (nimage > 0){ dipolesign = 1.0; }
		      else { dipolesign = -1.0; }
		      pz = phi->z[k] - (phi->zmin + 2.0 * (phi->zmax - phi->zmin) * double(nimage));
		      r3 = pow(r2 + pz*pz, 1.5);
		      deltaphi = p * dipolesign * pz / r3;
		      phi->data[index] += deltaphi;
		    }
		}
	    }
	}
    }
  return;
}

double MultiGrid::mu_Si (double E,double T)
{
  // Shamelssly copied from phosim
  // Jacobini et al. (1977) equation (9)
  double vm=1.53e9 * pow(T,-0.87); // cm/s
  double Ec = 1.01 * pow(T,1.55); // V/cm
  double beta = 2.57e-2 * pow(T,0.66); // index
  return((vm/Ec)/pow(1 + pow(fabs(E)/Ec,beta),1/beta));
}

void MultiGrid::FillRho(Array3D* rho, Array3D* elec, Array3D* hole)
{
  //Fill rho with data from electron array and hole array
  int i, j, k, index;
  double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (GridSpacing * GridSpacing);
  double ChargeDepth, TotalElectrons = 0.0, TotalHoles = 0.0;
  for (i=0; i<rho->nx; i++)
    { 
      for (j=0; j<rho->ny; j++)
	{
	  for (k=0; k<elec->nz; k++)
	    {
	      ChargeDepth = rho->Z(rho->zp[k] + rho->dzp / 2.0) - rho->Z(rho->zp[k] - rho->dzp / 2.0);  
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      rho->data[index] += (hole->data[index] - elec->data[index]) * RhoChargeFactor / ChargeDepth;
	      if (hole->data[index] < -1.0E-6 || elec->data[index] < -1.0E-6)
		{
		  printf("Negative hole or electron count! Something failed!\n");
		}
	      TotalElectrons += elec->data[index];
	      TotalHoles += hole->data[index];	      
	    }
	}
    }
  printf("Mobile charges added into rho.Total electrons=%.1f, Total holes=%.6g\n", TotalElectrons, TotalHoles);
  return;
}
