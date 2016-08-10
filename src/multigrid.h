/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** multigrid.h **************************

#include <stdio.h>            
#include <stdlib.h>           
#include <math.h>             
#include <time.h>
#include <string.h>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;

#include <globals.h>
#include <fileio.h>
#include <hdf5write.h>
#include <array3d.h>
#include <array2d.h>
#include <polygon.h>

class MultiGrid
{
 public:

  double w;		// Successive Over-Relaxation factor
  int ncycle;		// Number of SOR cycles at each resolution
  int iterations;		// Number of VCycles

  int ScaleFactor;       // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 1.0 micron, 96 grids in the z-direction
  double PixelSize;      // Pixel size in microns
  int GridsPerPixel;     // Number of grids per pixel at ScaleFactor = 1
  double GridSpacing;    // Grid size in microns
  int Nx;                // Number of grids in x at ScaleFactor = 1
  int Ny;                // Number of grids in y at ScaleFactor = 1
  int Nz;                // Number of grids in z at ScaleFactor = 1
  int Nzelec;                // Number of grids in z in electron arrays at ScaleFactor = 1  
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  double* SimulationRegionLowerLeft;	     
  int nsteps;            // Number of steps in Multigrid
  int XBCType;            // Free or periodic BC in X-direction
  int YBCType;            // Free or periodic BC in Y-direction  

  // Voltages and Charges
  double Vbb;		// Back bias
  double Vparallel_lo;	// Parallel Low Voltage
  double Vparallel_hi;	// Parallel High Voltage
  double Vserial_lo;	// Serial Low Voltage
  double Vaverage;       // Average voltage on bottom
  double Vchannelstop;   //  Voltage of undepleted channel stop
  double Vscupper;       // Scupper voltage

  int Channelkmin;              // Bottom of channel region doping 
  int Channelkmax;              // Top of channel region doping 
  int ChannelStopkmin;          // Bottom of channel stop region doping 
  int ChannelStopkmax;          // Top of channel stop region doping 
  int ChannelProfile;           // 0 = Square well, 1 = Gaussian
  int ChannelStopProfile;       // 0 = Square well, 1 = Gaussian
  double BackgroundDoping; 	// Background doping
  double ChannelStopDoping;	// Channel Stop doping
  double ChannelStopDepth;     	// Channel stop depth in microns
  double ChannelStopWidth;     	// Channel stop width in microns
  double ChannelStopCharge;
  double CSChargeDepth;
  double ChannelDoping;		// Channel doping
  double ChannelDepth;		// Channel depth in microns
  double GateOxide;             // Gate oxide thickness in microns
  int UndepletedChannelStop;	// 0 = No undepleted Region, 1 = undepleted Region

// Pixel Regions

  int NumberofPixelRegions;	  	  // 1
  double** PixelRegionLowerLeft;	  // 
  double** PixelRegionUpperRight;	  //
  int* NumberofFilledWells;		  //
  int** CollectedCharge;		  // Collected charge in e-
  double*** FilledPixelCoords;            // (x,y) coords of pixel center
  int CollectingPhases;                   // Number of collecting phases
  double CollectedChargeZmin;   	  // These parameters allow you to set the location of the initial collected charge
  double CollectedChargeZmax;
  double CollectedChargeXmin;
  double CollectedChargeXmax;
  double CollectedChargeYmin;
  double CollectedChargeYmax;
  
// Added dipoles

  int NumberofDipoles;                  // Number of dipoles
  double** DipoleCoords;                // X,Y coords of dipole center
  double* DipoleZLocation;               // Z location of dipole charge
  int* DipoleCharge;                     // Number of electrons in dipole
  int NumberofDipoleImages;             // Number of image pairs to calculate
  
// Constant Voltage Regions

  int NumberofFixedRegions;
  double** FixedRegionLowerLeft;
  double** FixedRegionUpperRight;
  double* FixedRegionVoltage;
  double* FixedRegionDoping;
  double* FixedRegionChargeDepth;
  int* FixedRegionBCType;  
  Array2D** BCType;      // BCType - 0->fixed potential; 1->Enormal = 0
  
  // Pixel Boundary Tests

  int PixelBoundaryTestType;
  int LogEField;
  int LogPixels;
  int LogPixelPaths;
  int PixelAreas;  
  double* PixelBoundaryLowerLeft;
  double* PixelBoundaryUpperRight;
  double* PixelBoundaryStepSize;

  int PixelBoundaryNx;
  int PixelBoundaryNy;
  int NumVertices;
  double ElectronZ0Area;
  double ElectronZ0Fill;
  double CCDTemperature;
  double DiffMultiplier;
  int NumDiffSteps;
  int SaturationModel;
  int NumElec;
  int NumSteps;
  double Sigmax;
  double Sigmay;
  double Xoffset;
  double Yoffset;    
  
  string outputfilebase; // Output filename base
  string outputfiledir; // Output filename directory

  int SaveData;
  int SaveElec;  
  
  Array3D** phi;      // Phi arrays
  Array3D** rho;      // Rho arrays
  Array3D** E;
  Array3D** elec;      // Number of stored electrons
  Array3D** hole;      // Number of mobile holes  
  
  MultiGrid() {};
  MultiGrid(string);
  ~MultiGrid();
 
  void ReadConfigurationFile(string);
  void BuildArrays(Array3D**, Array3D**, Array3D**, Array3D**, Array2D**);
  void SetInitialVoltages(Array3D*, Array2D*);
  void SetFixedCharges(Array3D*, Array2D*);
  void SetInitialElectrons(Array3D*, Array3D*);  
  void SetInitialHoles(Array3D*, Array3D*);
  void AdjustHoles(Array3D*, Array3D*, Array3D*);  
  void SOR(Array3D*, Array3D*, Array2D*, double);
  void SOR_N(Array3D*, Array3D*, Array2D*, double, int);  
  double Error(Array3D*, Array3D*);
  void Restrict(Array3D*, Array3D*, Array3D*, Array3D*, Array2D*, Array2D*);
  void Prolongate(Array3D*, Array3D*, Array2D*);
  void VCycle(Array3D**, Array3D**, Array2D**, double, int, int);
  void WriteOutputFile(string, string, string, Array3D*);
  void Gradient(Array3D*, Array3D**);
  void Trace(double*, int, bool, double, ofstream&);
  void TraceSpot(int);
  void TraceMultipleSpots(int);  
  void TraceGrid(int);  
  void TraceRegion(int);
  void FindEdge(double*, double, ofstream&);
  void FindCorner(double*, double*, ofstream&);
  void CalculatePixelAreas(int);
  void AddDipolePotentials(Array3D*);
  double mu_Si (double, double);
  void FillRho(Array3D*, Array3D*, Array3D*);  
};

