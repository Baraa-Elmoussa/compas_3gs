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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .egi import EGI
from .cell import Cell
from .polyhedral import FormNetwork
from .polyhedral import ForceVolMesh
from .disjointed import CellNetwork

__all__ = [
    "EGI",
    "Cell",
    "FormNetwork",
    "ForceVolMesh",
    "CellNetwork",
]
