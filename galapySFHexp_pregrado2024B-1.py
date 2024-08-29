#import emcee
#import dynesty
#import pytest
#import setuptools
#import requests
import numpy
import scipy
import matplotlib
import math, sys

import galapy

from galapy.analysis.plot import plt

#From galapy.internal.constants import clight


from galapy.Galaxy import GXY

ages=[1E+5,5E+5,1E+6,5E+6,1E+7,5E+7,1E+8,5E+8,1E+9,5E+9,1E+10] #yr,  time steps 

#choose galaxy ages
edad=1E+8
zz=0


#choose the time scale for an exponentially declining star formation history
tautau=1E+9

gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1,  csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+4, 'Zgxy' : 0.02},)


plt.title('exp SFH, tau=1E+9 yr, age=1E+8 yr')
xx=[edad, edad]
yy=[min(gxy.sfh(ages)), max(gxy.sfh(ages))]
plt.plot(xx, yy, 'b--') #to see where is the galaxy along the age vector

plt.plot(ages, gxy.sfh(ages), 'b-')

plt.xscale('log')
plt.xlabel('time [yr]')
plt.ylabel(r'SFR [M$_{\odot}$/yr]')

plt.show()


xx=[edad, edad]
yy=[min(gxy.sfh.Mstar(ages)), max(gxy.sfh.Mstar(ages))]
plt.plot(xx, yy, 'b--')

plt.plot(ages, gxy.sfh.Mstar(ages), 'b-')

print('the stellar mass accumulated till the galaxy age is', gxy.sfh.Mstar(edad), math.log10(gxy.sfh.Mstar(edad)))

plt.xscale('log')
plt.yscale('log')
plt.xlabel('time [yr]')
plt.ylabel(r'M* [M$_{\odot}$]')

plt.show()



waveEXP = gxy.wl() #array of rest-frame wavelengths in A
fluxnuEXP = gxy.get_SED() # Array of fluxes in milliJy (1E-3 microJy)

flambdaEXP=[]
for i in range(len(waveEXP)):
    flambdaEXP.append(fluxnuEXP[i]*1E+3*1E-29*(3E+18/waveEXP[i]**2))

plt.plot(waveEXP, flambdaEXP, 'b-')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel(r'f$_{\lambda}$ [erg/sec/cm^2/$\AA$]')

plt.show()





#compare SED of different ages

zz=0
tautau=1E+9


edad=1E+7
gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1,  csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+4, 'Zgxy' : 0.02},)


waveEXP7 = gxy.wl() #array of rest-frame wavelengths in A
fluxnuEXP7 = gxy.get_SED() # Array of fluxes in milliJy (1E-3 microJy)

flambdaEXP7=[]
for i in range(len(waveEXP7)):
    flambdaEXP7.append(fluxnuEXP7[i]*1E+3*1E-29*(3E+18/waveEXP7[i]**2))

plt.plot(waveEXP7, flambdaEXP7, 'b-', label='age=1E+7 yr')



edad=5E+8
gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1,  csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+4, 'Zgxy' : 0.02},)


waveEXP8 = gxy.wl() #array of rest-frame wavelengths in A
fluxnuEXP8 = gxy.get_SED() # Array of fluxes in milliJy (1E-3 microJy)

flambdaEXP8=[]
for i in range(len(waveEXP8)):
    flambdaEXP8.append(fluxnuEXP8[i]*1E+3*1E-29*(3E+18/waveEXP8[i]**2))

plt.plot(waveEXP8, flambdaEXP8, 'r-', label='age=5E+8 yr')

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel(r'f$_{\lambda}$ [erg/sec/cm^2/$\AA$]')



plt.show()





#compare SED for different tau



edad=1E+8


tautau=1E+9
gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1,  csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+4, 'Zgxy' : 0.02},)


waveEXP9 = gxy.wl() #array of rest-frame wavelengths in A
fluxnuEXP9 = gxy.get_SED() # Array of fluxes in milliJy (1E-3 microJy)

flambdaEXP9=[]
for i in range(len(waveEXP9)):
    flambdaEXP9.append(fluxnuEXP9[i]*1E+3*1E-29*(3E+18/waveEXP9[i]**2))


plt.title('SFR exp, age=1E+8 yr')
plt.plot(ages, gxy.sfh(ages), 'b-', label='tau=1E+9 yr')



tautau=1E+8
gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1,  csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+4, 'Zgxy' : 0.02},)


waveEXP8 = gxy.wl() #array of rest-frame wavelengths in A
fluxnuEXP8 = gxy.get_SED() # Array of fluxes in milliJy (1E-3 microJy)

flambdaEXP8=[]
for i in range(len(waveEXP8)):
    flambdaEXP8.append(fluxnuEXP8[i]*1E+3*1E-29*(3E+18/waveEXP8[i]**2))


plt.plot(ages, gxy.sfh(ages), 'r-', label='tau=1E+8 yr')

plt.legend()
plt.xscale('log')
plt.xlabel('time [yr]')
plt.ylabel(r'SFR [M$_{\odot}$/yr]')

plt.show()

    
plt.plot(waveEXP9, flambdaEXP9, 'b-', label='tau=1E+9 yr')
plt.plot(waveEXP8, flambdaEXP8, 'r-', label='tau=1E+8 yr')


plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel(r'f$_{\lambda}$ [erg/sec/cm^2/$\AA$]')


plt.show()





#make photometry

from galapy.PhotometricSystem import list_filters

#to choose all filters
#bands = list_filters('HST.WFC3.UVIS2')+list_filters('JWST.NIRcam')+list_filters('JWST.MIRI')+list_filters('ALMA')

#to select a few filters
pippo = ['HST.WFC3.UVIS2.F225W', 'HST.ACS.WFC.F625W']
bands = pippo+list_filters('JWST.MIRI')+list_filters('ALMA')
print(bands)

from galapy.PhotometricSystem import PMS
pms = PMS(*bands)
lpiv = pms.lpiv # wavelength of the filter configuration chosen


#to plot the filter transmission curves
from galapy.analysis.plot import photometric_system as pms_plot

fig, ax = plt.subplots(1,1,figsize=(12,3), constrained_layout=True)
_ = pms_plot(pms, ax=ax)
plt.show()


from galapy.Galaxy import PhotoGXY

edad=1E+8
tautau=1E+9
zz=0

gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1, csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+3, 'Zgxy' : 0.02},)

wave = gxy.wl() #array of rest-frame wavelengths in A
flux = gxy.get_SED() # array of fluxes in milliJy (1E-3 microJy)
waveZ=wave*(1+zz)
plt.plot(waveZ, flux, 'k-')


#to convolve with the filter transmission curves
pgxy = PhotoGXY(pms = pms, age = edad, redshift = zz, cosmo = 'Planck18', lstep=1, csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+3, 'Zgxy' : 0.02},)
pflux = pgxy.photoSED() #is this milliJ

print(pflux) #flux within a certain filter

plt.plot(lpiv, pflux, 'yD')



for i in range(len(lpiv)):
    print(lpiv[i], pflux[i], 23.9-2.5*math.log10(pflux[i]*1E+3)) 

plt.title('age=1E+8 yr, tau=1E+9 yr')
plt.ylabel(r'f${\nu}$ [m Jy]')
plt.xlabel(r'$\lambda$ [$\AA$]')

plt.xscale('log')
plt.yscale('log')
plt.show()




edad=1E+8
tautau=1E+9
zz=3

gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1, do_IGM = 'True', csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+3, 'Zgxy' : 0.02},)

wave = gxy.wl() #array of rest-frame wavelengths in A
flux = gxy.get_SED() # array of fluxes in milliJy (1E-3 microJy)
waveZ=wave*(1+zz) #to plot the spectrum in the observed frame
plt.plot(waveZ, flux, 'k-')


#to convolve with the filter transmission curves
pgxy = PhotoGXY(pms = pms, age = edad, redshift = zz, cosmo = 'Planck18', lstep=1, do_IGM = 'True', csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+3, 'Zgxy' : 0.02},)
pflux = pgxy.photoSED() #is this milliJ

print(pflux) #flux within a certain filter

plt.plot(lpiv, pflux, 'yD')



for i in range(len(lpiv)):
    print(lpiv[i], pflux[i], 23.9-2.5*math.log10(pflux[i]*1E+3)) 

plt.title('age=1E+8 yr, tau=1E+9 yr, z=3')
plt.ylabel(r'f${\nu}$ [m Jy]')
plt.xlabel(r'$\lambda$ [$\AA$]')

plt.xscale('log')
plt.yscale('log')
plt.show()




#to add the AGN contribution


edad=1E+8
tautau=1E+9
zz=3

gxy = GXY( age = edad, redshift = zz, cosmo = 'Planck18', lstep=1, do_IGM = 'True', do_AGN = 'True', csp = {'ssp_lib':'parsec22.NTL.refined'}, sfh = {'model': 'delayedexp', 'psi_norm' : 1., 'k_shape': 0, 'tau_star' : tautau, 'Mdust' : 1e+3, 'Zgxy' : 0.02},)


#To include Xray binaries we need   do_Xray = 'True'

wave = gxy.wl() #array of rest-frame wavelengths in A
flux = gxy.get_SED() # array of fluxes in milliJy (1E-3 microJy)
waveZ=wave*(1+zz) #to plot the spectrum in the observed frame
plt.plot(waveZ, flux, 'k-')

#here we define the components
components = gxy.components_to_flux()
print(list(components.keys()))
SEDstar=components['stellar']
SEDdust=components['extinct'] #stellar component extincted by dust without re emission 
SEDmc=components['MC'] #re emission of dust within molecular clouds
SEDdd=components['DD'] #diffuse dust re emission
#SEDxrb=components['XRB'] #xray binaries emission
SEDagn=components['AGN'] #contribution from AGN
plt.plot(wave, SEDstar, 'b-')
plt.plot(wave, SEDdust, 'g--')
plt.plot(wave, SEDmc, 'r-')
plt.plot(wave, SEDdd, 'm-')
#plt.plot(wave, SEDxrb, 'c-')
plt.plot(wave, SEDagn, 'y-')


plt.ylabel(r'f${\nu}$ [m Jy]')
plt.xlabel(r'$\lambda$ [$\AA$]')

plt.xscale('log')
plt.yscale('log')
plt.show()


sys.exit()
