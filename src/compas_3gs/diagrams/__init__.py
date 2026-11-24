"""
********************************************************************************
compas_3gs.diagrams
********************************************************************************

.. currentmodule:: compas_3gs.diagrams


Polyhedral Cell
===============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    EGI
    Cell

|

____


Multi-cell Polyhedron
=====================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ForceVolMesh
    FormNetwork

|

____


Cell Network
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    CellNetwork


"""

from .cell import Cell
from .disjointed import CellNetwork
from .egi import EGI
from .polyhedral import ForceVolMesh
from .polyhedral import FormNetwork

__all__ = [
    "EGI",
    "Cell",
    "FormNetwork",
    "ForceVolMesh",
    "CellNetwork",
]
