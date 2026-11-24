import pytest

from compas_3gs.algorithms import golden_section_search

# The search stops when its two interior points are closer than `tol`. The final interval is
# golden_ratio**3 = 2 + sqrt(5) times wider than that distance, it contains the minimum, and the
# result is its midpoint: the result is less than (2 + sqrt(5)) * tol / 2 from the minimum.
BOUND = (2 + 5**0.5) / 2

# the default `tol` of golden_section_search
TOL_DEFAULT = 1e-3


@pytest.mark.parametrize("xmin, a, b", [(2.0, 0.0, 5.0), (-1.3, -4.0, 3.0), (0.25, 0.0, 1.0)])
def test_golden_section_search_parabola(xmin, a, b):
    x = golden_section_search(lambda x: 3 * (x - xmin) ** 2 + 1, a, b)

    assert a < x < b
    assert abs(x - xmin) < BOUND * TOL_DEFAULT


@pytest.mark.parametrize("tol", [1e-3, 1e-6, 1e-9])
def test_golden_section_search_tol(tol):
    x = golden_section_search(lambda x: (x - 0.7) ** 2, 0.0, 1.0, tol=tol)

    assert abs(x - 0.7) < BOUND * tol


@pytest.mark.parametrize("function, xmin", [(lambda x: x, 1.0), (lambda x: -x, 3.0)])
def test_golden_section_search_minimum_at_an_end(function, xmin):
    x = golden_section_search(function, 1.0, 3.0)

    assert 1.0 < x < 3.0
    assert abs(x - xmin) < BOUND * TOL_DEFAULT
