import numpy as np

from fuvpy.src.background import _spencer_taper_weights


def test_spencer_taper_matches_reference_values():
    sza = np.array([0, 45, 90, 100, 120, 130, np.nan])
    taper = _spencer_taper_weights(sza)

    np.testing.assert_allclose(taper[:5], [1, 0.712256, 0.0685316, 0.00431824, 2.72559e-11])
    assert taper[5] == 0
    assert np.isnan(taper[6])


def test_spencer_taper_decreases_with_sza():
    taper = _spencer_taper_weights(np.arange(0, 121))

    assert np.all(np.diff(taper) <= 0)
