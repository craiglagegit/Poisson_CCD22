"""
build_spots.py

Modifying demo5.py to create 10 x 10 arrays of Gaussian spots
Craig Lage - 6-Apr-18

- Build a single large image, and access sub-images within it.
- Set the galaxy size based on the PSF size and a resolution factor.
- Set the object's flux according to a target S/N value.
- Shift by a random (dx, dy) drawn from a unit circle top hat.
"""

import sys, os, numpy, time
import logging
import galsim

def main(argv):
    """
    Make images to be used for characterizing the brighter-fatter effect
      - Each fits file is nx_tiles x ny_tiles postage stamps.
      - Each postage stamp is stamp_xsize x stamp_ysize pixels.
      - bf_nfile.fits, where nfile is just an identifier
    """
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("bf_plots")

    nx_tiles = 10                  #
    ny_tiles = 10                  #
    stamp_xsize = 20                #
    stamp_ysize = 20                #

    random_seed = 6424512           #

    pixel_scale = 0.2               # arcsec / pixel
    sky_level = 0.01                # ADU / arcsec^2

    gal_sigma = 0.3     # arcsec - sigma will be 1.5 pixels
    pixel_scale = 0.2  # arcsec / pixel

    shift_radius = 0.2  # arcsec

    logger.info('Starting bf_plots using:')
    logger.info('    - image with %d x %d postage stamps',nx_tiles,ny_tiles)
    logger.info('    - postage stamps of size %d x %d pixels',stamp_xsize,stamp_ysize)
    logger.info('    - Centroid shifts up to = %.2f pixels',shift_radius)

    sensor = galsim.Sensor()

    nfile = 1

    gal_file_name = 'bf_%d.fits'%nfile
    gal_flux = 2.0e5   # total counts on the image
    # Define the galaxy profile
    gal = galsim.Gaussian(flux=gal_flux, sigma=gal_sigma)
    logger.debug('Made galaxy profile')

    # This profile is placed at different locations
    # at each postage stamp in the gal image.
    gal_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)

    shift_radius_sq = shift_radius**2
    k=0
    starttime = time.time()

    for iy in range(ny_tiles):
        for ix in range(nx_tiles):
            ud = galsim.UniformDeviate(random_seed+k)

            b = galsim.BoundsI(ix*stamp_xsize+1 , (ix+1)*stamp_xsize-1, 
                               iy*stamp_ysize+1 , (iy+1)*stamp_ysize-1)
            sub_gal_image = gal_image[b]

            # Apply a random shift_radius:
            rsq = 2 * shift_radius_sq
            while (rsq > shift_radius_sq):
                dx = (2*ud()-1) * shift_radius
                dy = (2*ud()-1) * shift_radius
                rsq = dx**2 + dy**2

            this_gal = gal.shift(dx,dy)
            # Draw the image
            this_gal.drawImage(sub_gal_image)
            k += 1
    gal_image.write(gal_file_name)
    logger.info('Wrote image to %r',gal_file_name)  # using %r adds quotes around filename for us
    finishtime = time.time()
    print("Time to complete file = %.2f seconds\n"%(finishtime-starttime))
if __name__ == "__main__":
    main(sys.argv)
