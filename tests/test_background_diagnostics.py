import numpy as np

from fuvpy.src.background import backgroundmodel_BS
from test_background_dgsigma_freeze import _synthetic_wic


def test_fit_diagnostics_are_absent_by_default():
    result = backgroundmodel_BS(_synthetic_wic(), stop=1e-5)

    assert "dgfit_iterations" not in result
    assert "dgfit_converged" not in result


def test_fit_diagnostics_report_scalar_convergence():
    result = backgroundmodel_BS(
        _synthetic_wic(),
        stop=1e-5,
        diagnostics=True,
    )

    assert result["dgfit_iterations"].dims == ()
    assert result["dgfit_converged"].dims == ()
    assert int(result["dgfit_iterations"]) >= 2
    assert bool(result["dgfit_converged"])
    assert np.issubdtype(result["dgfit_iterations"].dtype, np.integer)
