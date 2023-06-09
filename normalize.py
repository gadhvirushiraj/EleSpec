# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 20:52:41 2022

@author: arvind
"""

import pandas as pd
import numpy as np
from scipy.ndimage import maximum_filter1d, median_filter

def spec_normalize(fl):
    '''
    Parameters
    ----------
    spec_df : Pandas DataFrame
        Pandas DataFrame containing the unnormalised spectrum. Must have columns 'wave' and 'flux'

    Returns
    -------
    norm_spec_df : Pandas DataFrame
        Pandas DataFrame containing the normalised spectrum.

    '''
    flux = np.array(fl)
    flux = flux/np.mean(flux) #Dividing all the flux values by the mean flux value for better accuracy.
    
    smooth_flux = flux
    for i in range(20): #Iterating to get a good continuum fit
        smooth_flux=maximum_filter1d(smooth_flux,int(30/(i+1)))
        smooth_flux=median_filter(smooth_flux,int(1000/(i+1))) #Required to conserve emission lines
    
    return flux/smooth_flux