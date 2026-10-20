#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 20:46:03 2022

@author: aohma
"""

import os
import warnings
import numpy as np
import pandas as pd
import xarray as xr

from datetime import datetime, timedelta

# from scipy.io import idl
from scipy.io import readsav

from fuvpy.src.measurement_uncertainty import calibration_multiplier, si_measurement_variance, wic_measurement_variance
from fuvpy.src.viewing_geometry import calculate_viewing_geometry
from fuvpy.utils.sunlight import subsol


DEFAULT_FRAME_QUALITY_THRESHOLDS = {
    'minimum_mapped_fraction': 0.0,
    'minimum_relative_spread': 0.1,
    'minimum_median_row_relative_spread': 0.06,
    'minimum_science_mlat': 68.0,
    'minimum_hv_mcp': 1000.0,
    'minimum_wic_hv_phos': 1000.0,
    'minimum_turn_on_mapped_fraction': 0.10,
    'maximum_turn_on_correlation': 0.90
}


def _viewing_state_valid(img):
    """Check the geometry fields needed to reconstruct this frame's rays."""
    position = np.asarray(img.O_GCI.values[0], dtype=float)
    latitude = float(img.LAT_CTR.values[0])
    longitude = float(img.LON_CTR.values[0])
    height = float(img.EMIS_HGT.values[0])
    return (np.all(np.isfinite(position)) and np.linalg.norm(position) > 6378.136
            and np.isfinite(latitude) and -90 <= latitude <= 90
            and np.isfinite(longitude) and np.isfinite(height) and height >= 0)


def _pointing_record_valid(img):
    """Flag categorical metadata errors without judging image or time-series quality."""
    axis = np.asarray(img.A_GCI3.values[0], dtype=float)
    spin = float(img.SPIN.values[0])
    phase = float(img.SPIN_PHASE.values[0])
    alignment = [float(img[name].values[0]) for name in ['INST_AZIMUTH', 'INST_CO_ELEV', 'INST_ROLL']]
    return bool(not np.isnat(img.pointing_time.values[0]) and _viewing_state_valid(img)
                and np.all(np.isfinite(axis)) and abs(np.linalg.norm(axis) - 1) <= 1e-3
                and np.isfinite(spin) and spin == np.rint(spin)
                and np.isfinite(phase) and np.all(np.isfinite(alignment)))


def read_idl(filenames, dzalim=80, hemisphere=None, reflat=True, remove_bad=True, measurement_variance=False, poisson_offset=0.0, wic_variance_model='reference', viewing_geometry=False, frame_quality_thresholds=None, pointing_metadata=False, skip_unmapped=False):
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
    viewing_geometry : bool, optional
        If True, reconstruct and store the pole-safe solar-view cosine and
        sine factors. The reconstructed DZA is checked against the DZA stored
        by FUVVIEW3. Default is False.
    pointing_metadata : bool, optional
        Retain spacecraft state, camera alignment, exact millisecond
        pointing_time and native coordinate metadata for shared WIC/SI12/SI13
        correction. Also enables viewing_geometry. Default False leaves
        existing output fields unchanged.
        Fit both background methods to a sensor mapping with
        pointing_correction=True to enable the shared repair.
    skip_unmapped : bool, optional
        During viewing-geometry validation, skip a file only when both stored
        and reconstructed DZA have no valid pixels. With pointing_metadata,
        retain such a frame as metadata-only when its pointing record is valid.
        Files with unreadable TIME are skipped with a warning in this mode.
        Default False preserves strict loading. Geometry disagreements raise.
    frame_quality_thresholds : dict, optional
        Overrides for ``minimum_mapped_fraction``,
        ``minimum_relative_spread``,
        ``minimum_median_row_relative_spread``, ``minimum_science_mlat``,
        detector high-voltage thresholds, and the maximum correlation allowed
        with a detector-off frame after turn-on. A turn-on correlation requires
        at least ``minimum_turn_on_mapped_fraction`` common mapped pixels.
        Quality 0 marks rejected frames, quality 1 marks detector-valid frames,
        and quality 2 marks science-ready frames.

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
    When ``viewing_geometry=True``:
    'solar_view_cosine_factor': Horizontal solar-viewing dot product
    'solar_view_sine_factor': Oriented horizontal solar-viewing cross product
    'sza': Solar zenith angle of each pixel
    'hemisphere': Hemisphere of the image
    'id': Short name of the detector
    'hv_phos': Phosphor high voltage stored by FUVVIEW3
    'hv_mcp': Microchannel-plate high voltage stored by FUVVIEW3
    'bad': Bool of problematic pixels. True=good and False=bad
    'frame_quality': Ordinal image classification: 0 rejected, 1 usable, and
                     2 science-ready
    When ``pointing_metadata=True``:
    'pointing_valid': Whether this frame has a usable pointing metadata record,
                      independently of image quality and mapped pixel count
    'row_relative_spread': Relative signal spread within each detector row
    'median_row_relative_spread': Median relative spread across mapped rows
    'turn_on_reference_correlation': Image correlation with the preceding
                                     detector-off frame during turn-on
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
    if pointing_metadata:
        viewing_geometry = True
    if isinstance(filenames, str):
        filenames = [filenames]

    filenames.sort()
    imgs = []
    skipped_unmapped = []
    skipped_invalid_time = []
    invalid_time_audit = []
    metadata_only = []
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
        try:
            image_time = pd.to_datetime(_timestamp_idl(imageinfo['time'][0]))
            if pointing_metadata:
                stamp = imageinfo['time'][0]
                milliseconds = float(stamp[1])
                if not np.isfinite(milliseconds) or milliseconds != int(milliseconds) or not 0 <= milliseconds < 86400000:
                    raise ValueError('invalid milliseconds in TIME')
                exact_time = datetime.strptime(str(int(stamp[0])), '%Y%j') + timedelta(milliseconds=int(milliseconds))
        except (KeyError, ValueError, TypeError, OverflowError) as error:
            if pointing_metadata and skip_unmapped:
                skipped_invalid_time.append(str(filenames[i]))
                invalid_time_audit.append(f'{filenames[i]}\tpointing_valid=false\treason=invalid_time\tdetail={error}')
                warnings.warn(f"Skipping {inst_id} frame {filenames[i]}: invalid TIME ({error})", RuntimeWarning, stacklevel=2)
                continue
            raise ValueError(f"Invalid TIME in {filenames[i]}: {error}") from error
        img = img.expand_dims(date=[image_time])
        img['hv_phos'] = (['date'], [float(imageinfo['hv_phos'][0])])
        img['hv_mcp'] = (['date'], [float(imageinfo['hv_mcp'][0])])

        if pointing_metadata:
            vectors = ['O_GCI', 'A_GCI3', 'SC_SV']
            scalars = ['SPIN_PHASE', 'SPIN', 'INST_AZIMUTH', 'INST_CO_ELEV', 'INST_ROLL', 'LAT_CTR', 'LON_CTR', 'EMIS_HGT']
            for name in vectors:
                try:
                    value = np.asarray(imageinfo[name.lower()][0], dtype=float).reshape(1, 3)
                except (KeyError, ValueError, TypeError, IndexError):
                    value = np.full((1, 3), np.nan)
                img[name] = (['date', 'pointing_xyz'], value)
            for name in scalars:
                try:
                    value = float(imageinfo[name.lower()][0])
                except (KeyError, ValueError, TypeError, IndexError):
                    value = np.nan
                img[name] = ('date', [value])
            img['pointing_time'] = ('date', [np.datetime64(exact_time, 'ns')])
            img['pointing_valid'] = ('date', [_pointing_record_valid(img)])
            for name in ['GLAT', 'GLON', 'SZA', 'DZA']:
                img['pointing_'+name] = (['date', 'row', 'col'], np.asarray(imageinfo[name.lower()][0])[None].copy())

        coordinates = ['mlat','mlon','mlt','glat','glon','sza','dza']
        fillvals = [-1e+31,-1e+31,-1e+31,-1e+31,-1e+31,-1,-1]

        # Replace fill values with np.nan
        for coordinate,fillval in zip(coordinates, fillvals):
            img[coordinate] = xr.where(img[coordinate]==fillval,np.nan,img[coordinate])

        # FUVVIEW3 stores the spacecraft radius in O_GCI and its Earth-fixed
        # direction in LAT_CTR and LON_CTR. Reproducing the stored DZA checks
        # that the reconstructed line of sight follows the same convention.
        if viewing_geometry:
            required_fields = ['o_gci', 'lat_ctr', 'lon_ctr', 'emis_hgt']
            available_fields = imageinfo.dtype.names if hasattr(imageinfo, 'dtype') and imageinfo.dtype.names else imageinfo.keys()
            available_fields = {field.lower() for field in available_fields}
            missing_fields = [field.upper() for field in required_fields if field not in available_fields]
            if missing_fields and not pointing_metadata:
                raise ValueError(f"Viewing geometry requires FUVVIEW3 field(s): {', '.join(missing_fields)}")

            geometry_ready = not pointing_metadata or _viewing_state_valid(img)
            if geometry_ready:
                spacecraft_radius = np.linalg.norm(np.asarray(imageinfo['o_gci'][0], dtype=float))
                spacecraft_latitude = float(imageinfo['lat_ctr'][0])
                spacecraft_longitude = float(imageinfo['lon_ctr'][0])
                emission_height = float(imageinfo['emis_hgt'][0])
                subsolar_latitude, subsolar_longitude = subsol([image_time])
                geometry = calculate_viewing_geometry(img['glat'].values, img['glon'].values, spacecraft_radius, spacecraft_latitude, spacecraft_longitude, emission_height, subsolar_latitude[0], subsolar_longitude[0])
            else:
                empty = np.full(img['dza'].shape, np.nan)
                geometry = {'dza': empty, 'solar_view_cosine_factor': empty.copy(), 'solar_view_sine_factor': empty.copy()}

            if geometry_ready:
                valid_dza = np.isfinite(geometry['dza']) & np.isfinite(img['dza'].values) & (img['dza'].values >= 0)
                if not np.any(valid_dza):
                    stored_valid = np.isfinite(img['dza'].values) & (img['dza'].values >= 0)
                    reconstructed_valid = np.isfinite(geometry['dza']) & (geometry['dza'] >= 0)
                    if skip_unmapped and not np.any(stored_valid) and not np.any(reconstructed_valid):
                        if pointing_metadata and bool(img.pointing_valid.values[0]):
                            metadata_only.append(str(filenames[i]))
                        else:
                            skipped_unmapped.append(str(filenames[i]))
                            warnings.warn(f"Skipping unmapped {inst_id} frame {filenames[i]}: no valid stored or reconstructed DZA pixels", RuntimeWarning, stacklevel=2)
                            continue
                    else:
                        raise ValueError("Viewing geometry could not be validated because the frame contains no valid DZA pixels")

                if np.any(valid_dza):
                    maximum_dza_error = np.max(np.abs(geometry['dza'][valid_dza] - img['dza'].values[valid_dza]))
                    if maximum_dza_error > 1e-3:
                        raise ValueError(f"Reconstructed viewing geometry disagrees with stored DZA by up to {maximum_dza_error:.4g} degrees")

            img['solar_view_cosine_factor'] = (['date', 'row', 'col'], geometry['solar_view_cosine_factor'])
            img['solar_view_sine_factor'] = (['date', 'row', 'col'], geometry['solar_view_sine_factor'])
            coordinates.extend(['solar_view_cosine_factor', 'solar_view_sine_factor'])

        # Select hemisphere
        mapped_mlat = img['mlat'].values[np.isfinite(img['mlat'].values)]
        median_mlat = np.median(mapped_mlat) if mapped_mlat.size else np.nan
        if hemisphere =='south' or (hemisphere is None and median_mlat<0):
            img=img.assign({'hemisphere': 'south'})
        elif hemisphere == 'north' or (hemisphere is None and median_mlat>0):
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

    if not imgs:
        raise ValueError("No frames remain after excluding unusable IDL files")
    imgs = xr.concat(imgs, dim='date', data_vars='all')
    imgs = imgs.assign({'id':  imageinfo['inst_id'][0].strip().decode('utf8')})
    if skipped_unmapped:
        imgs.attrs['skipped_unmapped_files'] = '\n'.join(skipped_unmapped)
    if skipped_invalid_time:
        imgs.attrs['skipped_invalid_time_files'] = '\n'.join(skipped_invalid_time)
        imgs.attrs['pointing_invalid_time_audit'] = '\n'.join(invalid_time_audit)
    if metadata_only:
        imgs.attrs['metadata_only_files'] = '\n'.join(metadata_only)
    imgs = _add_frame_quality(imgs, frame_quality_thresholds)
    if pointing_metadata:
        imgs.attrs['pointing_dzalim'] = float(dzalim)
        imgs.attrs['pointing_metadata'] = 'Original FUVVIEW3 state and native coordinates retained; correction not yet applied'
        imgs.O_GCI.attrs = {'long_name': 'Original spacecraft position in GCI', 'units': 'km'}
        imgs.pointing_time.attrs = {'long_name': 'Original IMAGE timestamp including milliseconds'}
        imgs.pointing_valid.attrs = {'long_name': 'Usable per-frame pointing metadata', 'definition': 'finite spacecraft state, unit spin axis, integer spin, valid camera geometry fields; independent of image quality and mapped pixels', 'spin_axis_norm_tolerance': 1e-3}

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
    imgs['hv_phos'].attrs = {'long_name': 'FUVVIEW3 phosphor high voltage', 'units': 'V'}
    imgs['hv_mcp'].attrs = {'long_name': 'FUVVIEW3 microchannel-plate high voltage', 'units': 'V'}
    if viewing_geometry:
        imgs['solar_view_cosine_factor'].attrs = {'long_name': 'Cosine-like solar-viewing azimuth factor', 'units': '1', 'definition': 'horizontal_solar dot horizontal_viewing'}
        imgs['solar_view_sine_factor'].attrs = {'long_name': 'Sine-like solar-viewing azimuth factor', 'units': '1', 'definition': 'local_up dot (horizontal_solar cross horizontal_viewing)'}

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
        'poisson_offset': float(poisson_offset),
        'calculated_from': 'img before optional WIC reflattening',
        'excludes': 'calibration uncertainty, background-model discrepancy, and atmospheric-response uncertainty'}
    if sensor == 'WIC':
        imgs['img_variance'].attrs['wic_variance_model'] = wic_variance_model
        imgs['img_variance'].attrs['reflattening_assumption'] = 'WIC reflattening is treated as a deterministic additive adjustment; img_variance is retained unchanged'

    return imgs


def _add_frame_quality(imgs, thresholds=None):
    """Classify frames from detector state, mapping, signal, and MLAT reach."""

    settings = dict(DEFAULT_FRAME_QUALITY_THRESHOLDS)
    if thresholds is not None:
        if not isinstance(thresholds, dict):
            raise TypeError("frame_quality_thresholds must be a dictionary or None")
        unknown = set(thresholds)-set(settings)
        if unknown:
            raise ValueError(f"Unknown frame-quality threshold(s): {', '.join(sorted(unknown))}")
        settings.update(thresholds)

    for name, value in settings.items():
        if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or not np.isfinite(value):
            raise ValueError(f"{name} must be a finite number")
    if not 0 <= settings['minimum_mapped_fraction'] <= 1:
        raise ValueError("minimum_mapped_fraction must lie between 0 and 1")
    if settings['minimum_relative_spread'] < 0 or settings['minimum_median_row_relative_spread'] < 0:
        raise ValueError("relative-spread thresholds must be non-negative")
    if not 0 <= settings['minimum_science_mlat'] <= 90:
        raise ValueError("minimum_science_mlat must lie between 0 and 90 degrees")
    if settings['minimum_hv_mcp'] < 0 or settings['minimum_wic_hv_phos'] < 0:
        raise ValueError("detector high-voltage thresholds must be non-negative")
    if not 0 <= settings['minimum_turn_on_mapped_fraction'] <= 1:
        raise ValueError("minimum_turn_on_mapped_fraction must lie between 0 and 1")
    if not 0 <= settings['maximum_turn_on_correlation'] <= 1:
        raise ValueError("maximum_turn_on_correlation must lie between 0 and 1")

    sensor = str(imgs['id'].values)
    image = imgs['img'].transpose('date', 'row', 'col').values
    detector_operational = np.isfinite(imgs['hv_mcp'].values) & (imgs['hv_mcp'].values >= settings['minimum_hv_mcp'])
    if sensor == 'WIC':
        detector_operational &= np.isfinite(imgs['hv_phos'].values) & (imgs['hv_phos'].values >= settings['minimum_wic_hv_phos'])

    mlat = imgs['mlat'].transpose('date', 'row', 'col').values
    sza = imgs['sza'].transpose('date', 'row', 'col').values
    dza = imgs['dza'].transpose('date', 'row', 'col').values
    mapped_pixels = np.isfinite(image) & np.isfinite(mlat) & np.isfinite(sza) & np.isfinite(dza)
    quality = np.zeros(imgs.sizes['date'], dtype=np.int8)
    row_relative_spread = np.full((imgs.sizes['date'], imgs.sizes['row']), np.nan)
    median_row_relative_spread = np.full(imgs.sizes['date'], np.nan)

    for frame in range(imgs.sizes['date']):
        mapped = mapped_pixels[frame]
        mapped_fraction = np.mean(mapped)
        if not np.any(mapped) or mapped_fraction < settings['minimum_mapped_fraction']:
            continue

        for row in range(imgs.sizes['row']):
            row_values = image[frame, row, mapped[row]]
            if row_values.size < 8:
                continue

            row_mean = np.abs(np.mean(row_values))
            row_spread = np.std(row_values)
            row_relative_spread[frame, row] = row_spread/row_mean if row_mean > np.finfo(float).eps else (0 if row_spread == 0 else np.inf)

        finite_rows = np.isfinite(row_relative_spread[frame])
        if not np.any(finite_rows):
            continue

        median_row_relative_spread[frame] = np.median(row_relative_spread[frame, finite_rows])

        if not detector_operational[frame]:
            continue

        values = image[frame][mapped]
        mean = np.abs(np.mean(values))
        spread = np.std(values)
        relative_spread = spread/mean if mean > np.finfo(float).eps else (0 if spread == 0 else np.inf)
        if not np.isfinite(relative_spread) or relative_spread < settings['minimum_relative_spread']:
            continue

        if median_row_relative_spread[frame] < settings['minimum_median_row_relative_spread']:
            continue

        quality[frame] = 1
        if np.max(np.abs(mlat[frame][mapped])) >= settings['minimum_science_mlat']:
            quality[frame] = 2

    # After a recorded detector off-to-on transition, reject operational frames
    # that still reproduce the detector-off pattern. Other quality failures do
    # not start a correlation chain.
    turn_on_reference_correlation = np.full(imgs.sizes['date'], np.nan)
    turn_on_frames = np.flatnonzero(detector_operational[1:] & ~detector_operational[:-1]) + 1
    for start in turn_on_frames:
        reference = start - 1
        if np.mean(mapped_pixels[reference]) < settings['minimum_turn_on_mapped_fraction']:
            continue

        for frame in range(start, imgs.sizes['date']):
            if not detector_operational[frame]:
                break

            common_mapping = mapped_pixels[reference] & mapped_pixels[frame]
            if np.mean(common_mapping) < settings['minimum_turn_on_mapped_fraction']:
                break

            correlation = _image_correlation(image[reference], image[frame], common_mapping)
            turn_on_reference_correlation[frame] = correlation
            if not np.isfinite(correlation) or correlation <= settings['maximum_turn_on_correlation']:
                break

            quality[frame] = 0

    imgs['frame_quality'] = (['date'], quality)
    imgs['row_relative_spread'] = (['date', 'row'], row_relative_spread)
    imgs['median_row_relative_spread'] = (['date'], median_row_relative_spread)
    imgs['turn_on_reference_correlation'] = (['date'], turn_on_reference_correlation)
    imgs['row_relative_spread'].attrs = {
        'long_name': 'Relative signal spread within each mapped detector row',
        'definition': 'standard deviation divided by absolute mean',
        'minimum_pixels': 8
    }
    imgs['median_row_relative_spread'].attrs = {
        'long_name': 'Median relative signal spread across mapped detector rows',
        'minimum_required': float(settings['minimum_median_row_relative_spread'])
    }
    imgs['turn_on_reference_correlation'].attrs = {
        'long_name': 'Correlation with the preceding detector-off frame during turn-on',
        'maximum_allowed': float(settings['maximum_turn_on_correlation']),
        'minimum_common_mapped_fraction': float(settings['minimum_turn_on_mapped_fraction']),
        'pixels': 'finite pixels mapped in both frames',
        'rejection_rule': 'after detector turn-on, reject frames while correlation exceeds threshold'
    }
    imgs['frame_quality'].attrs = {
        'long_name': 'Frame usability classification',
        'flag_values': np.array([0, 1, 2], dtype=np.int8),
        'flag_meanings': 'rejected usable science_ready',
        'minimum_mapped_fraction': float(settings['minimum_mapped_fraction']),
        'minimum_relative_spread': float(settings['minimum_relative_spread']),
        'minimum_median_row_relative_spread': float(settings['minimum_median_row_relative_spread']),
        'minimum_science_mlat': float(settings['minimum_science_mlat']),
        'minimum_hv_mcp': float(settings['minimum_hv_mcp']),
        'minimum_wic_hv_phos': float(settings['minimum_wic_hv_phos']),
        'minimum_turn_on_mapped_fraction': float(settings['minimum_turn_on_mapped_fraction']),
        'maximum_turn_on_correlation': float(settings['maximum_turn_on_correlation']),
        'detector_state_requirement': 'WIC requires HV_PHOS and HV_MCP; other sensors require HV_MCP',
        'science_mlat_definition': 'maximum absolute mapped magnetic latitude',
        'turn_on_rule': 'compare operational frames with the preceding detector-off frame after an observed off-to-on transition'
    }
    return imgs


def _image_correlation(reference, image, mask=None):
    """Pearson correlation between finite pixels in two detector images."""

    valid = np.isfinite(reference) & np.isfinite(image)
    if mask is not None:
        valid &= mask
    if np.count_nonzero(valid) < 2:
        return np.nan

    reference_values = reference[valid] - np.mean(reference[valid])
    image_values = image[valid] - np.mean(image[valid])
    denominator = np.sqrt(np.sum(reference_values**2)*np.sum(image_values**2))
    if denominator <= np.finfo(float).eps:
        return np.nan

    return np.sum(reference_values*image_values)/denominator

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
        The scalar background is estimated only from frames with
        ``frame_quality >= 1`` and the correction is applied to every frame.
    '''
    basepath = os.path.dirname(__file__)
    flatfields = readsav(basepath+'/../utils/wic_flatfield_dbase.idl')['flatfields']
    if pd.to_datetime(wic['date'][0].values)<pd.to_datetime('2000-10-03 23:30'):
        flat = flatfields[:,0]
    else:
        flat = flatfields[:,1]

    usable_frames = wic['frame_quality'].values >= 1
    background_pixels = (wic['sza'].values > 100) | (np.isnan(wic['sza'].values) & (wic['img'].values > 0))
    background_pixels &= usable_frames[:, None, None]
    background=np.nanmedian((wic[inImg]/flat[None,:,None]).values[background_pixels])
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
