#!/usr/bin/env python

import os
import os.path

import numpy as np
import numpy.lib.stride_tricks

import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm


def read_hdf5_array(f_name, verbose=True):
    with h5py.File(f_name, 'r') as f:
        if verbose:
            print 'Opened {0}'.format(f_name)
        # There should only be one item in this file.
        assert len(f.items()) == 1, 'Expected only one dataset.'
        ds_name, ds = f.items().pop()
        if verbose:
            print 'Loaded {0} with dimensions {1}.'.format(ds_name, ds.shape)
        return np.array(ds)


def block_view(A, block=(3, 3)):
    """Provide a 2D block view to 2D array."""
    assert len(A.shape) == 2 and len(block) == 2, 'Bad input shapes'
    assert A.shape[0] % block[0] == 0, 'Bad block[0]'
    assert A.shape[1] % block[1] == 0, 'Bad block[1]'
    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved ;-)
    shape= (A.shape[0] / block[0], A.shape[1] / block[1])+ block
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    return numpy.lib.stride_tricks.as_strided(A, shape= shape, strides= strides)


def pixel_sum(A):
    projected = np.sum(A, axis=2)
    blocked = block_view(projected[16:304, 16:304], block=(32, 32))
    return np.sum(blocked, axis=(-2, -1))


def load_grid(path):
    x_grid = np.loadtxt(os.path.join(path, 'grid_x.dat'), dtype=float)
    y_grid = np.loadtxt(os.path.join(path, 'grid_y.dat'), dtype=float)
    z_grid = np.loadtxt(os.path.join(path, 'grid_z.dat'), dtype=float)
    ze_grid = np.loadtxt(os.path.join(path, 'grid_ze.dat'), dtype=float)
    return x_grid, y_grid, z_grid, ze_grid


def read_step(
    index=20, band='y', name='BF_256_9x9_{band}_{index}', path=None,
    npix=3, pix_size=10., grid_size=110., nslab=4, save=True,
    verbose=False, log_color=True, V0=5.75):

    xy_lo, xy_hi = -0.5 * pix_size * npix, 0.5 * pix_size * npix

    # Load coordinate grids.
    grid_path = os.path.join(path, band)
    try:
        x_grid, y_grid, z_grid, ze_grid = load_grid(path=grid_path)
    except Exception as e:
        raise RuntimeError('Cannot read grids from {0}'.format(grid_path))
    x_center = 0.5 * (x_grid[1:] + x_grid[:-1]) - 0.5 * grid_size
    y_center = 0.5 * (y_grid[1:] + y_grid[:-1]) - 0.5 * grid_size
    z_center = 0.5 * (z_grid[1:] + z_grid[:-1])
    ze_center = 0.5 * (ze_grid[1:] + ze_grid[:-1])
    x_offset = x_grid - 0.5 * grid_size
    y_offset = y_grid - 0.5 * grid_size

    midx, midy = len(x_grid) // 2 - 1, len(y_grid) // 2 - 1
    assert np.allclose([x_grid[midx], y_grid[midy]], 0.5 * grid_size)
    slabx = slice(midx - nslab, midx + nslab)
    slaby = slice(midy - nslab, midy + nslab)
    if verbose:
        print 'Slab thickness: {0:.3f}um (x), {1:.3f}um (y).'.format(
            x_grid[midx + nslab] - x_grid[midx - nslab],
            y_grid[midy + nslab] - y_grid[midy - nslab])

    # Read outputs from this step into arrays.
    base_name = os.path.join(path, band, name.format(band=band, index=index))
    try:
        elec = read_hdf5_array(base_name + '_Elec.hdf5', verbose)
        hole = read_hdf5_array(base_name + '_Hole.hdf5', verbose)
        phi = read_hdf5_array(base_name + '_phi.hdf5', verbose)
    except IOError as e:
        raise RuntimeError('Cannot read {0}'.format(base_name))

    nelec, nhole = np.sum(elec), np.sum(hole)
    pix_elec = pixel_sum(elec)
    leakage = 100. * (1 - pix_elec[4, 4] / nelec)
    if verbose:
        print 'Mobile charges: {0} e + {1} h'.format(nelec, nhole)
        print 'Fraction outside central pixel: {0:.3f}%'.format(leakage)

    # Initialize the figure.
    fig = plt.figure(figsize=(14, 10.5))
    ax_z = plt.subplot2grid((2, 3), (0, 0))
    ax_xy = plt.subplot2grid((2, 3), (1, 0))
    ax_xz = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
    ax_yz = plt.subplot2grid((2, 3), (0, 2), rowspan=2)

    cmap = matplotlib.cm.get_cmap('viridis')
    cmap.set_over(cmap(1.))
    cmap.set_under(cmap(0.))
    cmap.set_bad(cmap(0.))
    cmap2 = matplotlib.cm.get_cmap('magma')
    cmap2.set_over(cmap2(1.))
    cmap2.set_under((0.,0.,0.,0.))
    cmap2.set_bad((0.,0.,0.,0.))

    pc_opts = dict(shading='flat', edgecolors='face')
    normalize = matplotlib.colors.LogNorm if log_color else matplotlib.colors.Normalize

    # === 1D plot along z
    phi_sum = np.mean(phi[slabx, slaby], axis=(0, 1))
    phi0_sum = np.mean(
        phi[midx + 64 - nslab: midx + 64 + nslab, midy + 64 - nslab: midy + 64 + nslab], axis=(0, 1))
    ax_z.plot(z_center, phi_sum, 'w-', lw=3)
    ax_z.plot(z_center, phi0_sum, 'w:', lw=3)

    t_opts = dict(fontsize=18, fontweight='bold', color='w',
                  horizontalalignment='right', verticalalignment='bottom')
    xyc = 'axes fraction'
    xy = (0.9, 0.3)
    label = '{0} band'.format(band)
    ax_z.annotate(label, xy, xy, xyc, xyc, **t_opts)
    xy = (0.9, 0.2)
    label = '{0:.1f}K e-'.format(1e-3 * nelec)
    ax_z.annotate(label, xy, xy, xyc, xyc, **t_opts)
    xy = (0.9, 0.1)
    label = 'leak {0:4.1f}%'.format(leakage)
    ax_z.annotate(label, xy, xy, xyc, xyc, **t_opts)

    ax_z_rhs = ax_z.twinx()
    elec_sum = np.sum(elec[slabx, slaby, :], axis=(0, 1))
    ax_z_rhs.hist(ze_center, bins=ze_grid, weights=elec_sum,
                  histtype='stepfilled', color=cmap(0.8), alpha=0.7)
    ax_z_rhs.set_ylim(0., None)
    ax_z_rhs.get_yaxis().set_ticks([])

    ax_z.set_axis_bgcolor(cmap(0.))
    ax_z.set_xlim(0., 2.5)
    ax_z.set_ylim(-1., 11.)
    ax_z.axhline(0., ls='--', c='w')
    ax_z.axhline(V0, ls='-', c='w')
    ax_z.set_xlabel('Distance from gates $z$ [um]')
    ax_z.set_ylabel('Central potential $\phi(z)$ [V]')

    # === 2D plot of bottom ~1.5um projected onto (x,y)
    for data, cm, alpha in ((elec, cmap, 1.0), (hole, cmap2, 0.5)):
        xy_proj = np.sum(data, axis=2)
        max_proj = np.percentile(xy_proj, 99.9)
        img = ax_xy.pcolormesh(
            x_offset, y_offset, xy_proj.T.copy(), cmap=cm, alpha=alpha,
            norm=normalize(vmax=max_proj, vmin=1e-3 * max_proj), **pc_opts)
    ax_xy.set_aspect(1.0)

    ax_xy.set_xlim(xy_lo, xy_hi)
    ax_xy.set_ylim(xy_lo, xy_hi)
    ax_xy.set_xlabel('Row (serial) offset $\Delta x$ [um]')
    ax_xy.set_ylabel('Column (parallel) offset $\Delta y$ [ym]')

    # === 2D plot of bottom ~1.5um in (x,z) slab.
    for data, cm, alpha in ((elec, cmap, 1.0), (hole, cmap2, 1.0)):
        xz_proj = np.sum(data[:, slaby, :], axis=1)
        max_proj = np.percentile(xz_proj, 99.9)
        img = ax_xz.pcolormesh(
            x_offset, ze_grid, xz_proj.transpose().copy(), cmap=cm, alpha=alpha,
            norm=normalize(vmax=max_proj, vmin=1e-3 * max_proj), **pc_opts)

    phi_sum = np.mean(phi[:, slaby, :], axis=1)
    ax_xz.contour(y_center, z_center, phi_sum.T.copy(),
                  levels=(0., V0), colors=('w', 'w'), linestyles=('--', '-'),
                  linewidths=(1, 1), alpha=1.0)

    ax_xz.set_xlim(xy_lo, xy_hi)
    ax_xz.set_ylim(ze_grid[0], ze_grid[-1])
    ax_xz.set_xlabel('Row (serial) offset $\Delta x$ [um]')
    ax_xz.set_ylabel('Distance from gates $z$ [um]')

    # === 2D plot of bottom ~1.5um in (y,z) slab.
    yz_elec = np.sum(elec[slabx, :, :], axis=0)
    max_elec = np.percentile(yz_elec, 99.9)
    img = ax_yz.pcolormesh(
        y_offset, ze_grid, yz_elec.T, cmap=cmap,
        norm=normalize(vmax=max_elec, vmin=1e-3 * max_elec), **pc_opts)

    phi_sum = np.mean(phi[slabx, :, :], axis=0)
    ax_yz.contour(y_center, z_center, phi_sum.T.copy(),
                  levels=(0., V0), colors=('w', 'w'), linestyles=('--', '-'),
                  linewidths=(1, 1), alpha=1.0)

    ax_yz.set_xlim(xy_lo, xy_hi)
    ax_yz.set_ylim(ze_grid[0], ze_grid[-1])
    ax_yz.set_xlabel('Column (parallel) offset $\Delta y$ [ym]')
    ax_yz.set_ylabel('Distance from gates $z$ [um]')

    opts = dict(c='w', alpha=0.75, ls=':')
    for ipix in range(1, npix):
        xy = xy_lo + ipix * pix_size
        ax_xy.axvline(xy, **opts)
        ax_xy.axhline(xy, **opts)
        ax_xz.axvline(xy, **opts)
        ax_yz.axvline(xy, **opts)

    plt.tight_layout()
    if save:
        plt.savefig('plots/{band}-{index}.png'.format(band=band, index=index))
    plt.close('all')


def main():

    path = os.environ['POISSON_BANDS']
    for band in 'ugrizy':
        for index in range(0, 110, 10):
            print('Processing {0} {1}'.format(band, index))
            try:
                read_step(band=band, index=index, path=path)
            except Exception as e:
                print(str(e))


if __name__ == '__main__':
    main()
