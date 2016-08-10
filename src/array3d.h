/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** array3d.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        

#include <globals.h>


class Array3D //This packages the 3D data sets
{
 public:
  int nx, ny, nz, nkz, *kz;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dzp, volume, *x, *y, *z, *zp, *zplus, *zminus, *dzpdz, *zpint, *data;
  Array3D() {};
  Array3D(double, double, int, double, double, int, double, double, int);
  ~Array3D();
  double PyramidalKernel3D(double, double, double);
  double DataInterpolate3D(double, double, double);
  int XIndex(double);
  int YIndex(double);
  int ZIndex(double);  
  double ZP(double);
  double ZPInt(double);
  double DZPDz(double);
  double D2ZPDz2(double);
  double Z(double);
  double ZPlus(double);
  double ZMinus(double);
};
