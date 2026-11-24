"""
********************************************************************************
compas_3gs.datastructures
********************************************************************************

.. currentmodule:: compas_3gs.datastructures


network3gs
==========

inherits compas mesh and adds these functions...

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Network3gs

|

____


mesh3gs
=======

inherits compas mesh and adds these functions...

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Mesh3gs

|

____


volmesh3gs
==========

inherits compas mesh and adds these functions...

.. autosummary::
    :toctree: generated/
    :nosignatures:

    VolMesh3gs

"""

from .mesh3gs import Mesh3gs
from .network3gs import Network3gs
from .volmesh3gs import VolMesh3gs

__all__ = [
    "Mesh3gs",
    "Network3gs",
    "VolMesh3gs",
]
