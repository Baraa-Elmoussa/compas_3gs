"""
********************************************************************************
compas_3gs.operations
********************************************************************************

.. currentmodule:: compas_3gs.operations


"""

from .celloperations import cell_split_indet_face_vertices
from .celloperations import cell_collapse_short_edge
from .celloperations import cell_relocate_face
from .celloperations import cell_face_subdivide_barycentric
from .celloperations import cell_merge_coplanar_adjacent_faces
from .volmeshoperations import point_reflection
from .volmeshoperations import volmesh_vertex_lift
from .volmeshoperations import volmesh_vertex_merge
from .volmeshoperations import volmesh_halfface_pinch
from .volmeshoperations import volmesh_halfface_extrude
from .volmeshoperations import volmesh_halfface_subdivide
from .volmeshoperations import volmesh_merge_adjacent_halffaces
from .volmeshoperations import volmesh_cell_split_vertex
from .volmeshoperations import volmesh_cell_subdivide_barycentric

__all__ = [
    "cell_split_indet_face_vertices",
    "cell_collapse_short_edge",
    "cell_relocate_face",
    "cell_face_subdivide_barycentric",
    "cell_merge_coplanar_adjacent_faces",
    "point_reflection",
    "volmesh_vertex_lift",
    "volmesh_vertex_merge",
    "volmesh_halfface_pinch",
    "volmesh_halfface_extrude",
    "volmesh_halfface_subdivide",
    "volmesh_merge_adjacent_halffaces",
    "volmesh_cell_split_vertex",
    "volmesh_cell_subdivide_barycentric",
]
