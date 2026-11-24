"""
********************************************************************************
compas_3gs.algorithms
********************************************************************************

.. currentmodule:: compas_3gs.algorithms


EGI
===

.. autosummary::
    :toctree: generated/
    :nosignatures:

    egi_from_vectors


____


Duality
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    volmesh_dual_volmesh
    volmesh_dual_network


____


Reciprocation
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    volmesh_reciprocate


____


Planarisation
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    volmesh_planarise


Arearisation
============


.. autosummary::
    :toctree: generated/
    :nosignatures:

    cell_arearise_face


____


Unified Diagram
===============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    volmesh_ud
    cellnetwork_ud

"""


from .arearisation import cell_arearise_face
from .duality import volmesh_dual_volmesh
from .duality import volmesh_dual_network
from .egi import egi_from_vectors
from .egi import cell_from_egi
from .other import golden_section_search
from .planarisation import cell_planarise
from .planarisation import volmesh_planarise
from .reciprocation import volmesh_reciprocate
from .unifieddiagram import volmesh_ud
from .unifieddiagram import cellnetwork_ud


__all__ = [
    "cell_arearise_face",
    "volmesh_dual_volmesh",
    "volmesh_dual_network",
    "egi_from_vectors",
    "cell_from_egi",
    "golden_section_search",
    "cell_planarise",
    "volmesh_planarise",
    "volmesh_reciprocate",
    "volmesh_ud",
    "cellnetwork_ud",
]
