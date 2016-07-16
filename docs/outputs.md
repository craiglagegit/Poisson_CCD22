# Output Files

This page documents the output files produced by the Poisson program and read by the various python scripts.

All output files are either in HDF5 format (with extension `.hdf5`) or plain text (with extension `.dat`).

## Charge density

File names ending with `_rho.hdf5` contain a 3D array of charge density values expressed `rho/eps` in units of Volts/(micron^2), where `eps` is the permittivity of silicon.

## Electric potential

File names ending with `_phi.hdf5` contain a 3D array of electric potential values in units of Volts.

## Electric field

File names ending with `_Ex.hdf5`, `_Ey.hdf5`, or `_Ez.hdf5` contain 3D arrays of electric field component (x,y,z) values in units of Volts/micron.

## Mobile charges

File names ending with `_Elec.hdf5` or `_Hole.hdf5` contain 3D arrays of mobile electron and hole counts in each grid cell. In other words, summing these arrays gives the total number of electrons and holes.

## Electron Paths

File names ending with `_Pts.dat` contain 3D coordinates of electrons traced through the sensor under the influence of the electric field and thermal diffusion.  The number of electrons in the file is set by the `NumElec` parameter and their initial distribution is determined by `PixelBoundaryTestType`.  When `LogPixelPaths = 0`, only the initial and final vertices are logged.  Otherwise, the vertices of each individual scatter are logged, leading to a much bigger file.
