"""
********************************************************************************
compas_3gs.utilities
********************************************************************************

.. currentmodule:: compas_3gs.utilities


Geometry
================================================================================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    resultant_vector
    datastructure_centroid


Topology
================================================================================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    pair_hf_to_uv
    pair_uv_to_hf


"""

from .geometry import cell_face_areaness
from .geometry import cell_face_flatness
from .geometry import datastructure_centroid
from .geometry import polygon_area_footprint
from .geometry import polygon_area_oriented
from .geometry import polygon_normal_oriented
from .geometry import resultant_vector
from .geometry import scale_polygon
from .geometry import volmesh_face_areaness
from .geometry import volmesh_face_flatness
from .other import compare_initial_current
from .other import get_force_colors_hf
from .other import get_force_colors_uv
from .other import get_force_mags
from .other import get_index_colordict
from .other import print_result
from .other import valuedict_to_colordict
from .topology import pair_hf_to_uv
from .topology import pair_uv_to_hf

__all__ = [
    "resultant_vector",
    "polygon_normal_oriented",
    "polygon_area_oriented",
    "polygon_area_footprint",
    "scale_polygon",
    "cell_face_flatness",
    "cell_face_areaness",
    "volmesh_face_flatness",
    "volmesh_face_areaness",
    "datastructure_centroid",
    "get_index_colordict",
    "valuedict_to_colordict",
    "compare_initial_current",
    "get_force_mags",
    "get_force_colors_uv",
    "get_force_colors_hf",
    "print_result",
    "pair_hf_to_uv",
    "pair_uv_to_hf",
]
