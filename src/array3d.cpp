/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** array3d.cpp **************************

#include "array3d.h"

Array3D::Array3D(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny, double Zmin, double Zmax, int Nz)
{
  int i, j, k;
  nkz = 10000;
  nx = Nx; ny = Ny; nz = Nz; 
  xmin = Xmin; xmax = Xmax; ymin = Ymin; ymax = Ymax; zmin = Zmin; zmax = Zmax;
  dx = (xmax - xmin) / (double) nx;
  dy = (ymax - ymin) / (double) ny;
  dzp = (zmax - zmin) / (double) nz;
  volume = dx * dy * dzp;
  x = new double[nx]; y = new double[ny]; z = new double[nz], zp = new double[nz], zplus = new double[nz], zminus = new double[nz]; dzpdz = new double[nz]; kz = new int[nkz]; zpint = new double[nkz];
  data = new double[nx * ny * nz];
  for (k=0; k<nz; k++)
    {
      zp[k] = zmin + dzp/2.0 + (double) k * dzp;
      z[k] = Z(zp[k]);
      dzpdz[k] = DZPDz(z[k]);
      zplus[k] = ZPlus(z[k]);
      zminus[k] = ZMinus(z[k]);
      //printf("zp = %f, z = %f,k = %d, kz = %d\n",zp[k], z[k],k,kz[k]);
    }
  for (k=0; k<nkz; k++)
    {
      kz[k] = ZIndex((double)k / 100.0);
      zpint[k] = ZP((double)k / 100.0);
    }
  for (j=0; j<ny; j++)
    {
      y[j] = ymin + dy/2.0 + (double) j * dy;
    }
  for (i=0; i<nx; i++)
    {
      x[i] = xmin + dx/2.0 + (double) i * dx;	  
      for (j=0; j<ny; j++)
	{
	  for (k=0; k<nz; k++)
	    {
	      data[i + j*nx + k*nx*ny] = 0.0;
	    }
	}
    }
}

Array3D::~Array3D()
{
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] zp;
  delete[] zplus;
  delete[] zminus;
  delete[] dzpdz;
  delete[] data;
  delete[] kz;
  delete[] zpint; 
}

double Array3D::PyramidalKernel3D(double deltax, double deltay, double deltaz)
{
  if (deltax>=1.0 || deltay>=1.0 || deltaz>=1.0)
    {
      return 0.0;
    }
  else
    {
      return (1.0-deltax)*(1.0-deltay)*(1.0-deltaz);
    }
}

double Array3D::DataInterpolate3D(double xin, double yin, double zin)
{
  int i, j, k, m, n, ml, nl, p, pl;
  double d, norm, PK, deltax, deltay, deltaz; 
  i=(int)((xin - xmin) / dx);
  j=(int)((yin - ymin) / dy);
  k=kz[max(0,min(nkz-1,(int)(zin * 100.0)))];
  d=0.0; norm = 0.0;
  for (m=i-1; m<i+2; m++)
    {
      ml = max(0,min(nx-1,m));
      deltax=fabs((xin-x[ml])/dx);
      for (n=j-1; n<j+2; n++)
	{
	  nl = max(0,min(ny-1,n));
	  deltay=fabs((yin-y[nl])/dy);
	  for (p=k-1; p<k+2; p++)
	    {
	      pl = max(0,min(nz-1,p));
	      deltaz=fabs((ZPInt(zin)-zp[pl])/dzp);
	      PK = PyramidalKernel3D(deltax,deltay,deltaz);
	      norm += PK;
	      d += PK * data[ml + nl * nx + pl * nx * ny];
	    }
	}
    }
  if (norm < 1.0E-15)
    {
      return 0.0;
    }
  else
    {
      return d/norm;
    }
}

double Array3D::ZP(double z)
{
  double n = 10.0;
  return - 100.0 * (n - 1.0) * pow(z / 100.0, (n + 1.0)/n) + n * z;
}

double Array3D::ZPInt(double z)
{
  int index = max(0,min(nkz-2,(int)(z * 100.0)));
  return zpint[index] + (z - (double)index / 100.0) * 100.0 * (zpint[index + 1] - zpint[index]);
}

double Array3D::DZPDz(double z)
{
  double n = 10.0;
  return - (n - 1.0) * (n + 1.0) / n * pow(z / 100.0, 1.0 / n) + n;
}

double Array3D::D2ZPDz2(double z)
{
  double n = 10.0;
  return - (n - 1.0) * (n + 1.0) / (n * n) * pow(z / 100.0, 1.0 / n - 1.0) / 100.0;
}

double Array3D::Z(double zp)

// Inverts ZP(z) using Newton's method
{
  int i = 0;
  double error = 1.0, lastroot = zp, newroot;
  while (error>1e-12)
    {
      newroot = lastroot - (ZP(lastroot) - zp) / DZPDz(lastroot);
      error = fabs((newroot - lastroot) / lastroot);
      lastroot=newroot;
      i=i+1;
      if (i > 100)
	{
	  printf("Iterations exceeded in Z(zprime). Quitting\n");
	  return newroot;
	}	  
    }
  return newroot;
}

double Array3D::ZPlus(double z)
{
  return dx * dx / (dzp * dzp) * (pow(DZPDz(z), 2.0) + dzp / 2.0 * D2ZPDz2(z)); 
}

double Array3D::ZMinus(double z)
{
  return dx * dx / (dzp * dzp) * (pow(DZPDz(z), 2.0) - dzp / 2.0 * D2ZPDz2(z));   
}


int Array3D::XIndex(double x)
{
  return (int)((x - xmin) / dx);
}

int Array3D::YIndex(double y)
{
  return (int)((y - ymin) / dy);
}

int Array3D::ZIndex(double z)
{
  return max(0,min(nz-1,(int)((ZP(z) - zmin) / dzp)));
}

