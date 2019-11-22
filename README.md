*** This is no longer the most up-to-date code.  ***
Users should move to the master branch at https://github.com/craiglagegit/Poisson_CCD, which is more recent and has a larger suite of example configuration files.

********************************************************************************************************

# Poisson_CCD22 - 'hole20' branch
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 13-Oct-17

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 360^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  Installation instructions and instructions on how to run the many examples are in the README.pdf file.


Changes:

This 'hole20' branch is a fairly major revision, with a number of new features.  In addition, it has been tested much more extensively, and many of the examples listed in the next section have been verifed with actual data taken on the STA3800 detector.  The Poisson_30Sep17.pdf presentation in the docs directory has a presentation-style summary of these results.  The following is a list of the changes:

* Three different methods for adding electrons to the CCD pixels have now been included, and these are selected with the ElectronMethod parameter.  A description and comparison of these three methods is given in slides 64-68 of the Poisson_30Sep17.pdf presentation in the docs directory.

* A ChannelStopSideDiff parameter has been added to allow independent control of side diffusion of the channel stops.  The default is that this is equal to the FieldOxideTaper parameter, which is how it ran in the past.

* Several parameters have been added to allow tree rings to be added to the simulation.  The treering example below will make their usage clear.

* An option has been added to simulate fringes projected onto the CCD in addition to the spots that were used before.  This is selected with PixelBoundaryTestType = 3, as in the fringerun example below.

* The option to save all of the coarser multi-grids used in the multi-grid solver has been added with the SaveMultiGrids parameter.  This is useful when testing convergence. The smallpixel and smallcap examples below illustrate the usage.

* The plotting routines have been modified fairly extensively, and all placed in a pysrc subdirectory.  Most of the subroutines used by these plotting routines have been collected in the pysrc/pysubs.py file.  This has eliminated a lot of duplication that was present in the past.

* A GateGap parameter has been added to study the impact of gaps between the gate electrodes.  This should be considered experimental at this time, and is not included in any of the examples below.

Updates - 13Oct17:

* Fixed a couple of minor things uncovered when simulating the E2V sensor.

* The Xoffset and Yoffset are now added to the initial fringe electron positions.

* Added code which deals with electrons which reach the top of the silicon (where the light comes in).  An electron which reaches here is reflected, with some probability of absorption, which is set by the parameter TopAbsorptionProb.  Default value is 0.0.

Acknowledgements:
This research was supported by DOE grant DE-SC0009999 and NSF/AURA grant N56981C.