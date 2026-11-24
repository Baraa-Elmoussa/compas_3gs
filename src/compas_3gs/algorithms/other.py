
from collections.abc import Callable


__all__ = ['golden_section_search']


def golden_section_search(function: Callable[[float], float], a: float, b: float, tol: float = 1e-3) -> float:
    """Golden section search algorithm for finding the minimum or maximum of a strictly unimodal function
    by successively narrowing the range of values inside which the extremum is known to exist.

    Parameters
    ----------
    function
        Evaluation function.
    a
        Lower bound of the initial search interval.
    b
        Upper bound of the initial search interval.
    tol
        Convergence tolerance, an interval width: the search stops when the two interior points
        are closer than `tol`.

    Returns
    -------
    float
        The location of the minimum, the midpoint of the final interval.

    Notes
    -----
    Adapted after the Wikipedia article *Golden section search*.

    References
    ----------
    - Wikipedia *Golden section search*. Available at: <https://en.wikipedia.org/wiki/Golden-section_search>
    - Kiefer, J. (1953). Sequential minimax search for a maximum. In *Proceedings of the American Mathematical Society* 4(3).

    """

    gr = (5 ** 0.5 + 1) / 2

    c = b - (b - a) / gr
    d = a + (b - a) / gr

    while abs(c - d) > tol:

        if function(c) < function(d):
            b = d

        else:
            a = c

        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (b + a) / 2


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == '__main__':
    pass
