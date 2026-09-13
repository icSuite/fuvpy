from pathlib import Path

import numpy as np

from fuvpy.src.load import read_idl
from fuvpy.src.viewing_geometry import calculate_viewing_geometry


def test_solar_viewing_factors_are_independent_of_polar_longitude():
    geometry = calculate_viewing_geometry([90, 90], [0, 180], 50000, 60, 20, 130, 20, 30)

    np.testing.assert_allclose(geometry['solar_view_cosine_factor'][0], geometry['solar_view_cosine_factor'][1], atol=1e-14)
    np.testing.assert_allclose(geometry['solar_view_sine_factor'][0], geometry['solar_view_sine_factor'][1], atol=1e-14)


def test_bundled_idl_geometry_reproduces_stored_dza():
    for name in ['wic_data', 's12_data', 's13_data']:
        directory = Path(__file__).parents[1]/f'data/{name}'
        filename = sorted(directory.glob('*.idl'))[0]
        image = read_idl(str(filename), reflat=False, viewing_geometry=True)
        valid = np.isfinite(image['dza'].values)

        assert np.any(valid)
        assert np.all(np.isfinite(image['solar_view_cosine_factor'].values[valid]))
        assert np.all(np.isfinite(image['solar_view_sine_factor'].values[valid]))
        assert image['solar_view_cosine_factor'].shape == image['img'].shape
        assert image['solar_view_sine_factor'].shape == image['img'].shape
