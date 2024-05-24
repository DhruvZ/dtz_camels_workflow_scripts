import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants



gnum = 78
nref_high = 32
nref_low = 1
base = '/orange/narayanan/d.zimmerman/camels_scripts/'
#========================================================
run = f'{base}snap33.galaxy{gnum}_nref{nref_high}.rtout.sed'
run2 = f'{base}snap33.galaxy{gnum}_nref{nref_low}.rtout.sed'
#========================================================


def sed_load(filename):
        m = ModelOutput(filename)
        wav,flux = m.get_sed(inclination='all',aperture=-1)
        wav  = np.asarray(wav)*u.micron #wav is in micron
        #wav *= (1.+z)
        flux = np.asarray(flux)[0]*u.erg/u.s
        return wav,flux

w1,f1 = sed_load(run)
w2,f2 = sed_load(run2)

print(f1)
print(f2)
print(f1/f2)
print(np.max((f1/f2)[(f1/f2).value>0]))
print(f1-f2)
plt.plot(w1,f1,label = f'nref{nref_high}',alpha=0.5)
plt.plot(w2,f2,label = f'nref{nref_low}',alpha=0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength')
plt.ylabel('flux')
plt.title(f'galaxy {gnum}')
plt.legend()
plt.savefig(f'/home/d.zimmerman/figures/s33_gal{gnum}_sedtest.png')
plt.close()


plt.plot(w1,f1/f2,label = f'nref{nref_high}/nref{nref_low}')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength')
plt.ylabel('flux ratio')
plt.title(f'galaxy {gnum}')
plt.legend()
plt.savefig(f'/home/d.zimmerman/figures/s33_gal{gnum}_sedtest_ratio.png')
plt.close()
