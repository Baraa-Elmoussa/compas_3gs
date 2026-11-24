from typing import Any
from typing import Self

from compas.datastructures import VolMesh
from compas.geometry import add_vectors
from compas_3gs.datastructures import Mesh3gs

__all__ = ["EGI"]


class EGI(Mesh3gs):
    """An EGI, represented as a mesh object.

    An EGI is a topological (mesh) dual of a cell of a volmesh, with all of its elements represented on a unit sphere.

    """

    def __init__(self, *args: Any, name: str = "EGI", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)

        self.attributes.update(
            {
                "origin": (0, 0, 0),
            }
        )
        self.default_vertex_attributes.update(
            {
                "x_fix": False,
                "y_fix": False,
                "z_fix": False,
                "type": None,
                "normal": None,
                "target_area": None,
            }
        )

    # --------------------------------------------------------------------------
    #   constructors
    # --------------------------------------------------------------------------

    @classmethod
    def from_volmesh_cell(cls, ckey: int, volmesh: VolMesh) -> Self:
        """Construct the EGI of a cell of a volmesh.

        Parameters
        ----------
        ckey
            The identifier of the cell.
        volmesh
            The volmesh, a [VolMesh][compas.datastructures.VolMesh].

        Returns
        -------
        Self
            The EGI, with a vertex at `origin + normal` for every halfface of the cell (vertex key = halfface key),
            and a face for every vertex of the cell (face key = vertex key), with outward normals.

        """
        egi = cls(name=str(ckey))
        origin = volmesh.cell_centroid(ckey)
        egi.attributes["origin"] = origin

        halffaces = volmesh.cell_faces(ckey)
        vertices = volmesh.cell_vertices(ckey)

        # mesh vertex j is vertices[j] and mesh face i is halffaces[i]
        mesh = volmesh.cell_to_mesh(ckey)

        for hfkey in halffaces:
            normal = volmesh.face_normal(hfkey)
            x, y, z = add_vectors(origin, normal)
            egi.add_vertex(key=hfkey, x=x, y=y, z=z)

        for j in mesh.vertices():
            face = [halffaces[f] for f in mesh.vertex_faces(j, ordered=True)][::-1]
            egi.add_face(face, fkey=vertices[j])

        return egi


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == "__main__":
    pass
