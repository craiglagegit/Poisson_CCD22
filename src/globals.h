/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** globals.h **************************

#define pi 4.0 * atan(1.0)     // Pi
#define max(x,y) (x>y?x:y)      // max macro definition
#define min(x,y) (x>y?y:x)      // min macro definition

#define QE             1.6E-19  // Electron charge in Coulombs
#define ME             9.11E-31 // Electron effective mass (approx = 1 in Si)
#define KBOLTZMANN     1.38E-23 // Boltzmann constant in Joules / Kelvin
#define EPSILON_0      8.85E-12 // Permittivity of free space in Farads / meter
#define EPSILON_SI     11.7     // Dielectric constant of silicon
#define EPSILON_OX     4.3      // Dielectric constant of silicon dioxide (should this be 3.9?)
#define MICRON_PER_M   1000000.0
#define MICRON_PER_CM  10000.0
#define METER_PER_CM   0.01
