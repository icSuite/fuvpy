#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 20:46:03 2022

@author: aohma
"""

import os
import numpy as np
import pandas as pd
import xarray as xr

from datetime import datetime

# from scipy.io import idl
from scipy.io import readsav

from fuvpy.src.measurement_uncertainty import calibration_multiplier, si_measurement_variance, wic_measurement_variance

def read_idl(filenames, dzalim=80, hemisphere=None, reflat=True, remove_bad=True, measurement_variance=False, poisson_offset=0.0, wic_variance_model='reference'):
    '''
    Load FUV images into a xarray.Dataset

    Parameters
    ----------
    filenames : str or list of str
        Path to one or several auroral fuv images stored in .idl or .sav files
    dzalim : float, optional
        Upper limit of the satellite zenith angle (viewing angle). The default is 80.
    hemisphere : str, optional
        Which hemisphere to include. 'north','south' and None.
        Default is None (automatically detected).
    reflat: bool, optional
        Reapply WIC's flatfield or keep original flatfield.
        Default is True (reapply).
    measurement_variance : bool, optional
        If True, add ``img_variance`` in the units of ``img`` squared. SI12
        and SI13 use reconstructed FUVVIEW3 calibration multipliers; WIC uses
        the published reference compound-Poisson approximation. Default is
        False, preserving the released loader behavior.
    poisson_offset : float, optional
        Nonnegative count added to the plug-in Poisson mean. Zero preserves
        the direct plug-in estimate. A value of 0.5 supplies an explicit
        low-count Jeffreys-Poisson offset. Used only when
        ``measurement_variance=True``.
    wic_variance_model : {'reference', 'correction_aware'}, optional
        WIC detector-variance model. ``reference`` uses the published 10.65 AD
        coefficient everywhere. ``correction_aware`` multiplies this
        coefficient by the reconstructed mission-gain and flat-field
        multiplier. Used only for WIC when ``measurement_variance=True``.
        Default is ``reference``.

    Returns
    -------
    xarray.Dataset with the following fields
    DATA:
    'img': Number of counts in each pixel
    'mlat': Magnetic latitude of each pixel
    'mlon': Magnetic longitude of each pixel
    'mlt': Magnetic local time of each pixel
    'glat': Geographic latitude of each pixel
    'glon': Geographic longitude of each pixel
    'dza': Satellite zenith angle (viewing angle) of each pixel
    'sza': Solar zenith angle of each pixel
    'hemisphere': Hemisphere of the image
    'id': Short name of the detector
    'bad': Bool of problematic pixels. True=good and False=bad
    When ``measurement_variance=True``:
    'img_variance': Conditional detector variance in img units squared
    'calibration_multiplier': FUVVIEW3 mission correction times flat field
    'calibration_flag': FUVVIEW3 calibration state for each image
    'airglow_scale': FUVVIEW3 airglow settings for each image
    DIMS:
    'date': Datetime of each image
    'row': rows of the detector
    'col': columns of the detector
    '''
    if isinstance(filenames, str):
        filenames = [filenames]

    filenames.sort()
    imgs = []
    calibration_flags = []
    airglow_scales = []
    for i in range(len(filenames)):
        imageinfo = readsav(filenames[i])['imageinfo']

        inst_id = imageinfo['inst_id'][0].strip().decode('utf8')

        # Read *.idl/*.sav files
        img = xr.Dataset({
        'img': (['row','col'],imageinfo['image'][0]),
        'mlat':(['row','col'],imageinfo['mlat'][0]),
        'mlon':(['row','col'],imageinfo['mlon'][0]),
        'mlt': (['row','col'],imageinfo['mlt'][0]),
        'glat':(['row','col'],imageinfo['glat'][0]),
        'glon':(['row','col'],imageinfo['glon'][0]),
        'dza': (['row','col'],imageinfo['dza'][0]),
        'sza': (['row','col'],imageinfo['sza'][0])
        })
        img = img.expand_dims(date=[pd.to_datetime(_timestamp_idl(imageinfo['time'][0]))])

        coordinates = ['mlat','mlon','mlt','glat','glon','sza','dza']
        fillvals = [-1e+31,-1e+31,-1e+31,-1e+31,-1e+31,-1,-1]

        # Replace fill values with np.nan
        for coordinate,fillval in zip(coordinates, fillvals):
            img[coordinate] = xr.where(img[coordinate]==fillval,np.nan,img[coordinate])

        # Select hemisphere
        if hemisphere =='south' or (hemisphere is None and np.nanmedian(img['mlat'])<0):
            img=img.assign({'hemisphere': 'south'})
        elif hemisphere == 'north' or (hemisphere is None and np.nanmedian(img['mlat'])>0):
            img=img.assign({'hemisphere': 'north'})
        else:
            img=img.assign({'hemisphere': 'na'})

        # coordinates to include
        ind = (img['dza'] < dzalim) & (img['mlat'] < 90) & (img['mlat'] > -90) & (img['glat']!=0)
        if inst_id=='WIC': ind = ind & (img['img'] > 0) # Zero-valued pixels in WIC are actually NaNs
        if remove_bad: ind = ind & ~_bad_pixels(inst_id) # Ignore bad parts of the detectors

        for coordinate in coordinates:
            img[coordinate] = xr.where(~ind,np.nan,img[coordinate])

        imgs.append(img)
        if measurement_variance:
            calibration_flags.append(int(imageinfo['calibration_flag'][0]))
            airglow_scales.append(np.asarray(imageinfo['airglow_scale'][0], dtype=float))

    imgs = xr.concat(imgs, dim='date')
    imgs = imgs.assign({'id':  imageinfo['inst_id'][0].strip().decode('utf8')})

    # Use the calibrated detector image before optional WIC reflattening.
    # Reflattening is treated as a deterministic additive adjustment, so this
    # conditional detector variance is retained unchanged afterward.
    if measurement_variance:
        imgs = _add_measurement_variance(imgs, calibration_flags, airglow_scales, poisson_offset, wic_variance_model)

    # Reapply WIC's flat field
    if (imgs['id']=='WIC')&reflat:
        imgs=_reflatWIC(imgs)

    # Add attributes
    imgs['img'].attrs = {'long_name': 'Original image'}
    imgs['mlat'].attrs = {'long_name': 'Magnetic latitude', 'units': 'deg'}
    imgs['mlon'].attrs = {'long_name': 'Magnetic longitude', 'units': 'deg'}
    imgs['mlt'].attrs = {'long_name': 'Magnetic local time', 'units': 'hrs'}
    imgs['glat'].attrs = {'long_name': 'Geographic latitude', 'units': 'deg'}
    imgs['glon'].attrs = {'long_name': 'Geographic longitude', 'units': 'deg'}
    imgs['sza'].attrs = {'long_name': 'Solar zenith angle', 'units': 'deg'}
    imgs['dza'].attrs = {'long_name': 'Viewing angle', 'units': 'deg'}

    return imgs


def _add_measurement_variance(imgs, calibration_flags, airglow_scales, poisson_offset, wic_variance_model='reference'):
    """Attach conditional detector variance to calibrated IMAGE/FUV data."""

    sensor = str(imgs['id'].values)
    calibration_flags = np.asarray(calibration_flags)
    airglow_scales = np.asarray(airglow_scales)

    if np.any(airglow_scales != 0):
        raise ValueError("Measurement variance is only implemented for IDL images without FUVVIEW3 airglow subtraction")

    if sensor in ['SI12', 'SI13']:
        if np.any(calibration_flags != 1):
            raise ValueError(f"{sensor} measurement variance requires calibration flag 1")
    elif sensor == 'WIC':
        # Files can carry flag 3 because the separate CIMAGE field is
        # background-subtracted. IMAGE is treated as the flag-1-equivalent
        # corrected-AD field used by fuvpy.
        if np.any(~np.isin(calibration_flags, [1, 3])):
            raise ValueError("WIC measurement variance requires calibration flag 1 or 3")
        if wic_variance_model not in ['reference', 'correction_aware']:
            raise ValueError("wic_variance_model must be 'reference' or 'correction_aware'")
    else:
        raise ValueError("Measurement variance is implemented for WIC, SI12, and SI13")

    multiplier = calibration_multiplier(imgs['date'].values, sensor)
    image = np.asarray(imgs['img'])

    if sensor in ['SI12', 'SI13']:
        variance = si_measurement_variance(image, multiplier, poisson_offset=poisson_offset)
        variance_model = ("SI Poisson counting variance: a * img + poisson_offset * a^2")
        variance_units = "corrected count squared"
    else:
        if wic_variance_model == 'reference':
            variance = wic_measurement_variance(image, poisson_offset=poisson_offset)
            variance_model = "WIC reference compound-Poisson variance: 10.65 * img + poisson_offset * 10.65^2"
        else:
            variance = wic_measurement_variance(image, poisson_offset=poisson_offset, multiplier=multiplier)
            variance_model = "WIC correction-aware compound-Poisson variance: (10.65 * a) * img + poisson_offset * (10.65 * a)^2"
        variance_units = "AD unit squared"

    valid = (np.isfinite(imgs['mlat']) & np.isfinite(multiplier) & (multiplier > 0))
    multiplier = np.where(valid, multiplier, np.nan)
    variance = np.where(valid, variance, np.nan)

    imgs['calibration_multiplier'] = (['date', 'row', 'col'], multiplier)
    imgs['img_variance'] = (['date', 'row', 'col'], variance)
    imgs['calibration_flag'] = (['date'], calibration_flags)
    imgs['airglow_scale'] = (['date', 'airglow_parameter'], airglow_scales)
    imgs = imgs.assign_coords(airglow_parameter=['multiplier', 'power', 'offset'])

    imgs['calibration_multiplier'].attrs = {'long_name': 'FUVVIEW3 mission correction times flat-field multiplier', 'units': '1'}
    imgs['img_variance'].attrs = {
        'long_name': 'Conditional detector measurement variance of img',
        'units': variance_units,
        'model': variance_model,
        'poisson_offset': float(poisson_offset)}
    if sensor == 'WIC':
        imgs['img_variance'].attrs['wic_variance_model'] = wic_variance_model

    return imgs

def _timestamp_idl(timestamp):
    """ returns datetime object for timestamp = imageinfo['TIME'][0] """

    hh = int(timestamp[1]/1000/60/60)
    mm = int(timestamp[1]/1000/60) - 60 * hh
    ss = int(timestamp[1]/1000) - 60**2 * hh - 60 * mm

    hh = str(hh).zfill(2)
    mm = str(mm).zfill(2)
    ss = str(ss).zfill(2)

    hhmmss = hh + mm + ss


    if len(str(timestamp[0])) != 7: # handle occasional errors in VIS image structure
        timestamp0 = int(str(timestamp[0])[:4] + str(timestamp[0])[4:].zfill(3))
    else:
        timestamp0 = timestamp[0]

    time = datetime.strptime(str(timestamp0) + hhmmss, '%Y%j%H%M%S')

    return time

def _bad_pixels(inst_id):
    '''
    Index of problematic pixels in the different detectors.
    The exact number will depent on the datetime, so this static approach is an approximation.
    Work could be done to get this function more accurate.

    Parameters
    ----------
    inst_id : str
        id of the FUV detector.

    Returns
    -------
    ind : np.array
        logical indices of bad pixels (True is bad).
    '''

    if inst_id == 'WIC':
        ind=np.zeros((256,256),dtype=bool)
        ind[230:,:] = True # Broken boom shades the upper rows of WIC (after what date?)
        ind[:2,:] = True #
    elif inst_id == 'SI13': #
        ind = np.zeros((128,128),dtype=bool)
        ind[:8,:] = True
        ind[120:,:] = True
    elif inst_id == 'VIS':
        ind = np.zeros((256,256),dtype=bool)
        ind[:,:4]=True
        # The corner of VIS Earth have an intensity dropoff
        for jj in range(25): ind[-1-jj,:25-jj] = True
        for jj in range(25): ind[-1-jj,-25+jj:] = True
        for jj in range(25): ind[jj,-25+jj:] = True
        for jj in range(25): ind[jj,:25-jj] = True
    elif inst_id=='UVI': # No known problems
        ind = np.zeros((228,200),dtype=bool)
    elif inst_id=='SI12': #
        ind = np.zeros((128,128),dtype=bool)
        ind[118:,:] = True # Broken boom shades the upper rows of WIC (after what date?)
        ind[:13,:] = True
    return ind

def _reflatWIC(wic,inImg='img',outImg='img'):
    '''
    The function changes how WIC's flatfield in implemented.
    Removes flatfield corrections, substractacts the general background noise
    from all the images, reapplies the flatfield and add background

    Parameters
    ----------
    imgs : xarray.Dataset
        Dataset with the FUV images.
    inImg : srtr, optional
        Name of the image to adjust. The default is 'image'.
    outImg : srtr, optional
        Name of the returned image. The default is 'image'.
    Returns
    -------
    imgs : xarray.Dataset
        Copy(?) of the FUV dataset with a new field containing the reflattened images.
    '''
    basepath = os.path.dirname(__file__)
    flatfields = readsav(basepath+'/../utils/wic_flatfield_dbase.idl')['flatfields']
    if pd.to_datetime(wic['date'][0].values)<pd.to_datetime('2000-10-03 23:30'):
        flat = flatfields[:,0]
    else:
        flat = flatfields[:,1]

    background=np.nanmedian((wic[inImg]/flat[None,:,None]).values[wic['sza'].values>100|(np.isnan(wic['sza'].values)&(wic['img'].values>0))])
    if np.isnan(background): background=450 # Set a reasonable value if no background pixels are available
    wic[outImg]=((wic[inImg].copy()/flat[None,:,None]-background)*flat[None,:,None]+background)
    return wic

def add_rayleigh(imgs,**kwargs):
    '''
    Add data_var with intinsity converted from counts to Rayleigh

    Parameters
    ----------
    imgs : xarray.Dataset
        Dataset with the FUV images.
    inImg : str, optional
        Name of the data_var with the input image in counts. The default is 'img'.
    inplace : bool, optional
        If True, update the Dataset in place. Default is False (return new Dataset) 
    
    Returns
    -------
    imgs : xarray.Dataset, optional
        Dataset including a new data_var with the intensity in kilo Rayleigh.
    '''
    
    # Set keyword arguments to input or default values    
    inImg = kwargs.pop('inImg') if 'inImg' in kwargs.keys() else 'img'
    inplace = bool(kwargs.pop('inplace')) if 'inplace' in kwargs.keys() else False

    # Make a copy if a new Dataset should be returned
    if not inplace: imgs = imgs.copy(deep=True)

    # Add a new data_var with counts in kilo Rayleigh
    if imgs['id']=='WIC':
        imgs[inImg+'R'] = imgs[inImg]/612.6
    elif imgs['id']=='VIS':
        imgs[inImg+'R'] = imgs[inImg]/4.32
    elif imgs['id']=='UVI':
        imgs[inImg+'R'] = imgs[inImg]/35.52
    elif imgs['id']=='SI13':
        imgs[inImg+'R'] = imgs[inImg]/15.3

    # Add attributes
    imgs[inImg+'R'].attrs = {'long_name': imgs[inImg].attrs['long_name'], 'units': 'kR'}
    
    # Return the new Dataset
    if not inplace: return imgs
