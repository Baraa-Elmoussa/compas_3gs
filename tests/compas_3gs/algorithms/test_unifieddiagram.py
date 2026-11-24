import json
import os

import pytest

from compas.datastructures import Mesh
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_ud
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork


@pytest.fixture
def forcediagram():
    """A 3x3x3 meshgrid of unit cubes."""
    return ForceVolMesh.from_meshgrid(dx=3, nx=3)


@pytest.fixture
def formdiagram(forcediagram):
    """The dual network of the meshgrid: one node at every cell centroid."""
    return volmesh_dual_network(forcediagram, cls=FormNetwork)


def test_volmesh_ud_keys(forcediagram, formdiagram):
    cells, prisms = volmesh_ud(forcediagram, formdiagram)

    assert sorted(cells) == sorted(forcediagram.cells())
    assert sorted(prisms) == sorted(formdiagram.edges())
    assert len(cells) == 27
    assert len(prisms) == 54


@pytest.mark.parametrize("scale", [0.25, 0.5, 0.75])
def test_volmesh_ud_cells(forcediagram, formdiagram, scale):
    cells, _ = volmesh_ud(forcediagram, formdiagram, scale=scale)

    for cell, scaled in cells.items():
        mesh = Mesh.from_vertices_and_faces(scaled["vertices"], scaled["faces"])
        assert mesh.number_of_vertices() == 8
        assert mesh.number_of_faces() == 6
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) == pytest.approx(scale**3)
        # every form node is at its cell centroid, so the scaled cell keeps the centroid
        assert TOL.is_allclose(centroid_points(scaled["vertices"]), forcediagram.cell_centroid(cell))


@pytest.mark.parametrize("scale", [0.25, 0.5, 0.75])
def test_volmesh_ud_prisms(forcediagram, formdiagram, scale):
    cells, prisms = volmesh_ud(forcediagram, formdiagram, scale=scale)

    for (u, v), prism in prisms.items():
        mesh = Mesh.from_vertices_and_faces(prism["vertices"], prism["faces"])
        assert len(prism["vertices"]) == 8
        assert sorted({index for face in prism["faces"] for index in face}) == list(range(8))
        assert mesh.is_closed()
        # a box between the scaled shared faces: scale x scale, and 1 - scale long
        assert volume_polyhedron(mesh.to_vertices_and_faces()) == pytest.approx(scale**2 * (1 - scale))
        # the first face is on the scaled cell u, the second on the scaled cell v
        for xyz in prism["vertices"][:4]:
            assert min(distance_point_point(xyz, other) for other in cells[u]["vertices"]) < TOL.absolute
        for xyz in prism["vertices"][4:]:
            assert min(distance_point_point(xyz, other) for other in cells[v]["vertices"]) < TOL.absolute


def test_volmesh_ud_moved_formdiagram(forcediagram, formdiagram):
    cells, prisms = volmesh_ud(forcediagram, formdiagram)
    for node in formdiagram.nodes():
        x, y, z = formdiagram.node_coordinates(node)
        formdiagram.node_attributes(node, "xyz", [x + 10, y - 5, z + 2])

    moved_cells, moved_prisms = volmesh_ud(forcediagram, formdiagram)

    # the form diagram is moved onto the force diagram first
    for cell in cells:
        assert moved_cells[cell]["faces"] == cells[cell]["faces"]
        for a, b in zip(cells[cell]["vertices"], moved_cells[cell]["vertices"]):
            assert distance_point_point(a, b) < TOL.absolute
    for edge in prisms:
        for a, b in zip(prisms[edge]["vertices"], moved_prisms[edge]["vertices"]):
            assert distance_point_point(a, b) < TOL.absolute


def test_volmesh_ud_honeycomb():
    with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
        data = json.load(f)
    forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])
    formdiagram = volmesh_dual_network(forcediagram, cls=FormNetwork)

    cells, prisms = volmesh_ud(forcediagram, formdiagram, scale=0.5)

    assert sorted(cells) == sorted(forcediagram.cells())
    for cell, scaled in cells.items():
        mesh = Mesh.from_vertices_and_faces(scaled["vertices"], scaled["faces"])
        original = forcediagram.cell_to_mesh(cell)
        assert mesh.number_of_vertices() == original.number_of_vertices()
        assert mesh.number_of_faces() == original.number_of_faces()
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) == pytest.approx(0.125 * volume_polyhedron(original.to_vertices_and_faces()))

    assert sorted(prisms) == sorted(formdiagram.edges())
    for prism in prisms.values():
        mesh = Mesh.from_vertices_and_faces(prism["vertices"], prism["faces"])
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0


@pytest.mark.parametrize("scale", [0, 1])
def test_volmesh_ud_scale_0_or_1_raises(forcediagram, formdiagram, scale):
    with pytest.raises(Exception, match="equivalent to the polyhedral"):
        volmesh_ud(forcediagram, formdiagram, scale=scale)


@pytest.mark.parametrize("scale", [-0.5, 1.5])
def test_volmesh_ud_scale_outside_raises(forcediagram, formdiagram, scale):
    with pytest.raises(AssertionError):
        volmesh_ud(forcediagram, formdiagram, scale=scale)
