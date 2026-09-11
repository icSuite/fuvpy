"""Reconstruct the IMAGE/FUV viewing direction needed by the BS model."""

import numpy as np


EARTH_EQUATORIAL_RADIUS_KM = 6378.136
EARTH_INVERSE_FLATTENING = 298.257


def calculate_viewing_geometry(glat, glon, spacecraft_radius, spacecraft_latitude, spacecraft_longitude, emission_height, subsolar_latitude, subsolar_longitude):
    """Calculate DZA and the solar-view factors for one IMAGE/FUV frame.

    The inputs are the geographic pixel coordinates and the FUVVIEW3 fields
    ``|O_GCI|``, ``LAT_CTR``, ``LON_CTR``, and ``EMIS_HGT``. Pixel positions
    are placed on the same emission ellipsoid used by FUVVIEW3.

    The two solar-view factors are the dot and oriented cross products of the
    horizontal solar and viewing vectors. Unlike a geographic azimuth, these
    definitions remain meaningful at the geographic pole.
    """

    glat = np.asarray(glat, dtype=float)
    glon = np.asarray(glon, dtype=float)
    latitude = np.deg2rad(glat)
    longitude = np.deg2rad(glon)

    # FUVVIEW3 increases both ellipsoid axes by the assumed emission height.
    equatorial_radius = EARTH_EQUATORIAL_RADIUS_KM + float(emission_height)
    polar_radius = EARTH_EQUATORIAL_RADIUS_KM*(1 - 1/EARTH_INVERSE_FLATTENING) + float(emission_height)
    eccentricity_squared = 1 - (polar_radius/equatorial_radius)**2
    prime_vertical_radius = equatorial_radius/np.sqrt(1 - eccentricity_squared*np.sin(latitude)**2)

    pixel = np.stack([
        prime_vertical_radius*np.cos(latitude)*np.cos(longitude),
        prime_vertical_radius*np.cos(latitude)*np.sin(longitude),
        prime_vertical_radius*(1 - eccentricity_squared)*np.sin(latitude)
    ], axis=-1)

    spacecraft_latitude = np.deg2rad(float(spacecraft_latitude))
    spacecraft_longitude = np.deg2rad(float(spacecraft_longitude))
    spacecraft = float(spacecraft_radius)*np.array([
        np.cos(spacecraft_latitude)*np.cos(spacecraft_longitude),
        np.cos(spacecraft_latitude)*np.sin(spacecraft_longitude),
        np.sin(spacecraft_latitude)
    ])

    line_of_sight = spacecraft - pixel
    line_of_sight /= np.linalg.norm(line_of_sight, axis=-1)[..., None]

    # FUVVIEW3 defines DZA relative to the geocentric pixel radius.
    up = pixel/np.linalg.norm(pixel, axis=-1)[..., None]
    cosine_dza = np.sum(line_of_sight*up, axis=-1)
    dza = np.rad2deg(np.arccos(np.clip(cosine_dza, -1, 1)))

    subsolar_latitude = np.deg2rad(float(subsolar_latitude))
    subsolar_longitude = np.deg2rad(float(subsolar_longitude))
    solar_direction = np.array([
        np.cos(subsolar_latitude)*np.cos(subsolar_longitude),
        np.cos(subsolar_latitude)*np.sin(subsolar_longitude),
        np.sin(subsolar_latitude)
    ])

    solar_horizontal = solar_direction - np.sum(solar_direction*up, axis=-1)[..., None]*up
    viewing_horizontal = line_of_sight - cosine_dza[..., None]*up

    return {
        'dza': dza,
        'solar_view_cosine_factor': np.sum(solar_horizontal*viewing_horizontal, axis=-1),
        'solar_view_sine_factor': np.sum(up*np.cross(solar_horizontal, viewing_horizontal), axis=-1)
    }
