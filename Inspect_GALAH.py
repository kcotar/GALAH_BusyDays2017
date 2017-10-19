
# coding: utf-8

# # Inspect_GALAH: A tool to for GALAH spectra as well as SME/Cannon output
# 
# Author: Sven Buder

# In[108]:

# Compatibility with Python 3
from __future__ import (absolute_import, division, print_function)

try:
    get_ipython().magic(u'matplotlib inline')
    get_ipython().magic(u"config InlineBackend.figure_format='retina'")
except:
    pass

# Basic packages
import numpy as np
import os
import sys

# Packages to work with FITS and (IDL) SME.out files
import astropy.io.fits as pyfits
import astropy.table as table
from scipy.io.idl import readsav

# Matplotlib and associated packages for plotting
import matplotlib.pyplot as plt
params = {'text.usetex': True, 'text.latex.preamble': [r'\usepackage{upgreek}', r'\usepackage{amsmath}'],'font.family' : 'lmodern','font.size' : 11}   
plt.rcParams.update(params)


# First be sure you are in the right directory and where your spectra and SME/Cannon files are:

# In[109]:

print('Your working directory is '+os.getcwd())


# Now we can start with the major routines to save the spectra and SME/Cannon files in a dictionary strucuture: 

# # Routine to get SME data in structure 'sme'

# In[110]:

def get_sme_data(directory='OUTPUT',field='bmstar',sobject_id=140708006401203,setup='DR2',mode='all',show_na=False):
    """
    INPUT:
    
    directory  = directory where the file is
    field      = field that was designated by WG4 when running SME
    sobject_id = sobject_id
    setup      = DR2 by default used by WG4 to initialise SME segments & line masks
    mode       = Instructions in 'mode_DR2' part of https://github.com/svenbuder/GALAH_BusyDays2017/blob/master/README.md
    show_na    = Show which modes could not be loaded, if mode=='all'
    
    OUTPUT:
    
    'sme'
    
    This is a dictionary structure, i.e. you can access the wavelength entry via sme['wave']. Keywords are:
    
    field, sobject_id, setup, mode : same as input

    Depending on what you use for 'mode', different arrays will be created with variable 'mode':    
    *mode*_wave : SME wavelength grid of given *mode*
    *mode*_sob  : SME norm. spectrum of given *mode* on 'wave' grid
    *mode*_uob  : SME error of sob of given *mode* on 'wave' grid
    *mode*_smod : SME's synthetic spectrum of given *mode* on 'wave' grid

    """
    
    sme = {}
    sme['field'] = field
    sme['sobject_id'] = sobject_id
    sme['setup'] = setup
    sme['mode'] = mode

    if mode == 'all':

        mode = [
            'Sp',
            'Li', 'C',
            'O', 'O7772', 'O7774', 'O7775', 'Na', 'Na4752', 'Na5683', 'Na5688',
            'Mg', 'Mg4730', 'Mg5711', 'Mg7692','Al', 'Al6696', 'Al6699', 'Al7835', 'Al7836',
            'Si', 'Si5666', 'Si5690', 'Si5793', 'Si6722', 'K', 'K5802',
            'Ca', 'Ca5857', 'Ca5862', 'Ca5868', 'Ca6494', 'Ca6500', 'Ca6501', 'Ca6509',
            'Sc', 'Ti', 'Fe', 'V', 'Cr', 'Mn', 'Co', 'Ni','Cu', 'Cu5700', 'Cu5782',
            'Zn', 'Zn4722', 'Zn4811', 'Rb', 'Sr', 'Y', 'Y4855', 'Y4884', 'Y5663', 'Y5729',
            'Zr', 'Mo', 'Ru', 'Ba', 'Ba5854', 'Ba6497', 'La', 'Ce', 'Nd', 'Nd4811',
            'Sm', 'Eu', 'Eu5819', 'Eu6645'
            ]

    for each_mode in mode:
        try:
            sme_out = readsav(directory+'/'+field+'_'+str(sobject_id)+'_'+setup+'_'+each_mode+'_SME.out')
            sme_out = sme_out.sme[0]

            sme[each_mode+'_wave'] = sme_out.wave
            sme[each_mode+'_sob'] = sme_out.sob
            sme[each_mode+'_uob'] = sme_out.uob
            sme[each_mode+'_smod'] = sme_out.smod
        except:
            if show_na | (len(mode)==1):
                print('Mode '+each_mode+' not available')
            else:
                pass
        
    return sme;


# Try this to call the function with default values:

# In[111]:

sme = get_sme_data()


# # Create 'fits' class and fill with necessary data and enable plot routines

# In[197]:

def get_fits_data(directory='default',sobject_id=140708006401203):
    """
    INPUT:
    
    directory  = directory where the 4 CCD fits files are
    sobject_id = sobject_id
    
    OUTPUT:
    
    'spectrum'
    
    This is a dictionary structure, i.e. you can access the wavelength entry via sme['wave']. Keywords are:
    
    directory  : same as input
    sobject_id : same as sobject_id
    wave_raw   : IRAF wavelength of raw flux
    sob_raw    : IRAF raw flux for CCD
    uob_raw    : IRAF raw flux error
    wave_norm  : IRAF wavelength of GUESS normalised flux
    sob_norm   : IRAF GUESS normalised flux
    uob_norm   : IRAF GUESS normalised flux error
    as well as ccd_*CCD* in front of the major keywords for CCD Nr. *CCD* (e.g. ccd_*CCD*_wave_raw)

    """
    
    spectrum = {}
    spectrum['directory'] = directory
    spectrum['sobject_id'] = sobject_id

    for each in range(1,4+1):

        if directory == 'default':
            door = pyfits.open('SPECTRA/dr5.2/'+str(sobject_id)[0:6]+'/standard/com/'+str(sobject_id)+str(each)+'.fits')
        else:
            door = pyfits.open(self.directory+str(sobject_id)+str(each)+'.fits')

        ws=door[0].header["CRVAL1"]
        inc=door[0].header["CDELT1"]
        nax=door[0].header["NAXIS1"]
        spectrum['ccd'+str(each)+'_wave_raw']=map(lambda x:((x+1)*inc+ws),range(0,nax))
        ws=door[4].header["CRVAL1"]
        inc=door[4].header["CDELT1"]
        nax=door[4].header["NAXIS1"]
        spectrum['ccd'+str(each)+'_wave_norm']=map(lambda x:((x+1)*inc+ws),range(0,nax))

        spectrum['ccd'+str(each)+'_sob_raw'] = door[0].data
        spectrum['ccd'+str(each)+'_sob_norm'] = door[4].data

        spectrum['ccd'+str(each)+'_uob_raw'] = door[0].data*door[1].data
        spectrum['ccd'+str(each)+'_uob_norm'] = door[4].data*door[1].data[len(door[1].data)-len(door[4].data):]

    spectrum['wave_raw'] = np.concatenate([spectrum['ccd'+str(ccd)+'_wave_raw'] for ccd in range(1,4+1)])
    spectrum['wave_norm'] = np.concatenate([spectrum['ccd'+str(ccd)+'_wave_norm'] for ccd in range(1,4+1)])
    spectrum['sob_raw'] = np.concatenate([spectrum['ccd'+str(ccd)+'_sob_raw'] for ccd in range(1,4+1)])
    spectrum['sob_norm'] = np.concatenate([spectrum['ccd'+str(ccd)+'_sob_norm'] for ccd in range(1,4+1)])
    spectrum['uob_raw'] = np.concatenate([spectrum['ccd'+str(ccd)+'_uob_raw'] for ccd in range(1,4+1)])
    spectrum['uob_norm'] = np.concatenate([spectrum['ccd'+str(ccd)+'_uob_norm'] for ccd in range(1,4+1)])

    return spectrum

class fits_class(object):
    """ A class with all FITS entries and necessary routines """
    
    def __init__(self, sobject_id=140708006401203, directory='default'):
        self.sobject_id = sobject_id
        self.directory = directory
        
        spectrum = get_fits_data(self.directory,self.sobject_id)
        
        self.wave_raw  = spectrum['wave_raw']
        self.ccd1_wave_raw  = spectrum['ccd1_wave_raw']
        self.ccd2_wave_raw  = spectrum['ccd2_wave_raw']
        self.ccd3_wave_raw  = spectrum['ccd3_wave_raw']
        self.ccd4_wave_raw  = spectrum['ccd4_wave_raw']
        self.wave_norm = spectrum['wave_norm']
        self.ccd1_wave_norm = spectrum['ccd1_wave_norm']
        self.ccd2_wave_norm = spectrum['ccd2_wave_norm']
        self.ccd3_wave_norm = spectrum['ccd3_wave_norm']
        self.ccd4_wave_norm = spectrum['ccd4_wave_norm']
        self.sob_raw   = spectrum['sob_raw']
        self.ccd1_sob_raw   = spectrum['ccd1_sob_raw']
        self.ccd2_sob_raw   = spectrum['ccd2_sob_raw']
        self.ccd3_sob_raw   = spectrum['ccd3_sob_raw']
        self.ccd4_sob_raw   = spectrum['ccd4_sob_raw']
        self.sob_norm  = spectrum['sob_norm']
        self.ccd1_sob_norm  = spectrum['ccd1_sob_norm']
        self.ccd2_sob_norm  = spectrum['ccd2_sob_norm']
        self.ccd3_sob_norm  = spectrum['ccd3_sob_norm']
        self.ccd4_sob_norm  = spectrum['ccd4_sob_norm']
        self.uob_raw   = spectrum['uob_raw']
        self.ccd1_uob_raw   = spectrum['ccd1_uob_raw']
        self.ccd2_uob_raw   = spectrum['ccd2_uob_raw']
        self.ccd3_uob_raw   = spectrum['ccd3_uob_raw']
        self.ccd4_uob_raw   = spectrum['ccd4_uob_raw']
        self.uob_norm  = spectrum['uob_norm']
        self.ccd1_uob_norm  = spectrum['ccd1_uob_norm']
        self.ccd2_uob_norm  = spectrum['ccd2_uob_norm']
        self.ccd3_uob_norm  = spectrum['ccd3_uob_norm']
        self.ccd4_uob_norm  = spectrum['ccd4_uob_norm']
        
    def plot_norm_spectrum_on4axes(self,
        include_errors=True,
        ylim=(-0.1,1.1),
        xlim_ccd4='raw',
        savefig='PLOTS_4CCDs',
        help_lines=False
        ):
        """
        INPUT:
        
        include_errors : Includes grey area sob-uob to sob+uob in background
        ylim           : (lower_yaxis,upper_yacis), default: (-0.1,1.1)
        savefig        : will save plot in directory passed with 'savefig'
        help_lines     : will add vertical lines with (name,lambda)
        
        OUTPUT:
        
        4 axes plot of normalised spectrum
        
        """
        
        kwargs_plot = dict(color='black', linewidth=0.5)
        kwargs_fill_between = dict(color='grey', linewidth=0.5)
        f,(ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(15,10),sharey=True)
        sob1 = ax1.plot(self.ccd1_wave_norm,self.ccd1_sob_norm,label='Spectrum',**kwargs_plot)
        ax2.plot(self.ccd2_wave_norm,self.ccd2_sob_norm,**kwargs_plot)
        ax3.plot(self.ccd3_wave_norm,self.ccd3_sob_norm,**kwargs_plot)
        ax4.plot(self.ccd4_wave_norm,self.ccd4_sob_norm,**kwargs_plot)
        if include_errors==True:
            fill1 = ax1.fill_between(self.ccd1_wave_norm,self.ccd1_sob_norm-self.ccd1_uob_norm,self.ccd1_sob_norm+self.ccd1_uob_norm,label='Error',**kwargs_fill_between)
            ax2.fill_between(self.ccd2_wave_norm,self.ccd2_sob_norm-self.ccd2_uob_norm,self.ccd2_sob_norm+self.ccd2_uob_norm,**kwargs_fill_between)
            ax3.fill_between(self.ccd3_wave_norm,self.ccd3_sob_norm-self.ccd3_uob_norm,self.ccd3_sob_norm+self.ccd3_uob_norm,**kwargs_fill_between)
            ax4.fill_between(self.ccd4_wave_norm,self.ccd4_sob_norm-self.ccd4_uob_norm,self.ccd4_sob_norm+self.ccd4_uob_norm,**kwargs_fill_between)
        ax1.legend(loc='lower left')
        ax4.set_xlabel(r'$\mathrm{Wavelength~[\AA]}$')
        ax1.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax2.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax3.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax4.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax1.set_ylim(ylim)
        if xlim_ccd4=='raw':
            ax4.set_xlim(self.ccd4_wave_raw[0],self.ccd4_wave_raw[-1])
        plt.tight_layout()
        
        if savefig!=False:
            plt.savefig(savefig+'/'+str(self.sobject_id)+'.pdf')
            
    def plot_raw_spectrum_on4axes(self,
        include_errors=True,
        ylim=(-0.1,1.1),
        savefig='PLOTS_4CCDs',
        help_lines=False
        ):
        """
        INPUT:
        
        include_errors : Includes grey area sob-uob to sob+uob in background
        ylim           : (lower_yaxis,upper_yacis), default: (-0.1,1.1)
        savefig        : will save plot in directory passed with 'savefig'
        help_lines     : will add vertical lines at given wavelengths
        
        OUTPUT:
        
        4 axes plot of raw spectrum
        
        """
        
        kwargs_plot = dict(color='black', linewidth=0.5)
        kwargs_fill_between = dict(color='grey', linewidth=0.5)
        f,(ax1,ax2,ax3,ax4) = plt.subplots(4,figsize=(15,10),sharey=True)
        sob1 = ax1.plot(self.ccd1_wave_norm,self.ccd1_sob_norm,label='Spectrum',**kwargs_plot)
        ax2.plot(self.ccd2_wave_norm,self.ccd2_sob_norm,**kwargs_plot)
        ax3.plot(self.ccd3_wave_norm,self.ccd3_sob_norm,**kwargs_plot)
        ax4.plot(self.ccd4_wave_norm,self.ccd4_sob_norm,**kwargs_plot)
        if include_errors==True:
            fill1 = ax1.fill_between(self.ccd1_wave_norm,self.ccd1_sob_norm-self.ccd1_uob_norm,self.ccd1_sob_norm+self.ccd1_uob_norm,label='Error',**kwargs_fill_between)
            ax2.fill_between(self.ccd2_wave_norm,self.ccd2_sob_norm-self.ccd2_uob_norm,self.ccd2_sob_norm+self.ccd2_uob_norm,**kwargs_fill_between)
            ax3.fill_between(self.ccd3_wave_norm,self.ccd3_sob_norm-self.ccd3_uob_norm,self.ccd3_sob_norm+self.ccd3_uob_norm,**kwargs_fill_between)
            ax4.fill_between(self.ccd4_wave_norm,self.ccd4_sob_norm-self.ccd4_uob_norm,self.ccd4_sob_norm+self.ccd4_uob_norm,**kwargs_fill_between)
        ax1.legend(loc='lower left')
        ax4.set_xlabel(r'$\mathrm{Wavelength~[\AA]}$')
        ax1.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax2.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax3.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax4.set_ylabel(r'$\mathrm{Flux~[norm]}$')
        ax4.set_xlim(s)
        
        if not help_lines:
            print('--')
        
        plt.tight_layout()
        
        if savefig!=False:
            plt.savefig(savefig+'/'+str(self.sobject_id)+'.pdf')


# # Now let's create the class FITS for a given sobject_id and use the provided functions on it!

# In[193]:

fits = fits_class(sobject_id=140708006401203)


# ## Assuming you want to look at specific lines (and check RV), then use define the 'help_lines' variable before plotting

# In[195]:

help_lines=[(r'$\mathrm{H_\alpha}$',6562.7970),(r'$\mathrm{H_\beta}$',4861.3230),(r'$\mathrm{Li}$',6707.7635)]


# ## Plot normalised spectrum on 4 axes

# In[196]:

fits.plot_norm_spectrum_on4axes(help_lines=help_lines)


# ## Plot raw spectrum on 4 axes

# In[169]:

fits.plot_raw_spectrum_on4axes()


# In[202]:

os.system('ipython nbconvert --to script Inspect_GALAH.ipynb');

