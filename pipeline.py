import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import os
from astropy.io import fits

DATADIR = '/home/andrej/projects/tess/data/sector01/'

def read_from_fits(tess_id, normalize=True):
    filename = glob.glob('%s/*%d*fits' % (DATADIR, tess_id))[0]

    with fits.open(filename) as hdu:
        # pylint: disable=E1101
        time = hdu[1].data['TIME']
        flux = hdu[1].data['PDCSAP_FLUX']
        flux_err = hdu[1].data['PDCSAP_FLUX_ERR']

        cond = np.where(np.logical_and(~np.isnan(time), ~np.isnan(flux), ~np.isnan(flux_err)))
        time = time[cond]
        if normalize:
            flux = flux[cond]/flux[cond].mean()
            flux_err = flux_err[cond]/flux.mean()
        else:
            flux = flux[cond]
            flux_err = flux_err[cond]

    return (time, flux, flux_err)

def bjd2phase(time, bjd0, period):
    phase = -0.5+((time-bjd0-period/2) % period) / period
    return phase

def run_lombscargle(time, flux, ferr, pmin=0.1, pmax=15, subsample=0.001, npeaks=3, extras=''):
    with tempfile.NamedTemporaryFile() as lcf:
        np.savetxt(lcf, np.vstack((time, flux, ferr)).T)

        cmdline = 'vartools -i %s -ascii -LS %f %f %f %d 1 /%s/ %s' % (lcf.name, pmin, pmax, subsample, npeaks, tempfile.gettempprefix(), extras)
        ls = subprocess.run(cmdline.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lsres = ls.stdout.split()

        freq, ls, logfap = np.loadtxt('%s.ls' % lcf.name, unpack=True)
        os.unlink('%s.ls' % lcf.name)
        
        rv = {}
        rv['freq'] = freq
        rv['ls'] = ls
        rv['logfap'] = logfap
        rv['periods'] = [float(lsres[1+4*i+0]) for i in range(npeaks)]
        rv['logfaps'] = [float(lsres[1+4*i+1]) for i in range(npeaks)]
        rv['persnrs'] = [float(lsres[1+4*i+2]) for i in range(npeaks)]
        rv['lsstats'] = [float(lsres[1+4*i+3]) for i in range(npeaks)]

        return rv

def process_lc(tess_id, bjd0=None, period=None):
    time, flux, flux_err = read_from_fits(tess_id)

    # light curve:
    plt.xlabel('Time')
    plt.ylabel('Normalized flux')
    plt.plot(time, flux, 'b.')
    plt.show()

    # run period finder(s):
    lsres = run_lombscargle(time, flux, flux_err, subsample=0.001)

    # plot periodogram(s):
    plt.xlabel('Period')
    plt.ylabel('LS')
    plt.plot(1./lsres['freq'], lsres['ls'], 'b-')
    for i in range(3):
        plt.axvline(lsres['periods'][i], ls='--', c='r')
    if period is not None:
        plt.axvline(period, ls='-', c='g')
    plt.show()

    plt.xlabel('Phase')
    plt.ylabel('Normalized flux')
    phase = bjd2phase(time, bjd0=bjd0, period=lsres['periods'][0])
    plt.plot(phase, flux, 'b.')
    plt.show()

# process_lc(tess_id=234518605, bjd0=1326.539962, period=5.675008)
process_lc(tess_id=55524055, bjd0=1325.391669, period=2.64942)