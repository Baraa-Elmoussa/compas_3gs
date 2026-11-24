import random

import pytest

from compas.geometry import angle_vectors
from compas.geometry import distance_point_point
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork

# the angle tolerances of volmesh_reciprocate, in degrees
TOLERANCE = 0.1
TOLERANCE_BOUNDARY = 1.0


@pytest.fixture
def forcediagram():
    """A 3x3x3 meshgrid of unit cubes."""
    return ForceVolMesh.from_meshgrid(dx=3, nx=3)


@pytest.fixture
def formdiagram(forcediagram):
    """The dual network of the meshgrid: every form edge is perpendicular to its face."""
    return volmesh_dual_network(forcediagram, cls=FormNetwork)


@pytest.mark.parametrize("weight", [0.0, 0.5, 1.0])
def test_volmesh_reciprocate_reciprocal_pair_stays(forcediagram, formdiagram, weight):
    force_xyz = {vertex: forcediagram.vertex_coordinates(vertex) for vertex in forcediagram.vertices()}
    form_xyz = {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()}

    volmesh_reciprocate(forcediagram, formdiagram, kmax=10, weight=weight, tolerance=TOLERANCE)

    assert max(distance_point_point(forcediagram.vertex_coordinates(vertex), xyz) for vertex, xyz in force_xyz.items()) < TOL.absolute
    assert max(distance_point_point(formdiagram.node_coordinates(node), xyz) for node, xyz in form_xyz.items()) < TOL.absolute


def test_volmesh_reciprocate_weight_1_perturbed_form(forcediagram, formdiagram):
    rng = random.Random(1)
    for node in sorted(formdiagram.nodes()):
        x, y, z = formdiagram.node_coordinates(node)
        formdiagram.node_attributes(node, "xyz", [x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), z + rng.uniform(-0.5, 0.5)])
    force_xyz = {vertex: forcediagram.vertex_coordinates(vertex) for vertex in forcediagram.vertices()}

    before = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )
    volmesh_reciprocate(forcediagram, formdiagram, kmax=50, weight=1.0, tolerance=TOLERANCE)
    after = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )

    assert before > 45
    assert after < TOLERANCE
    # the form edges keep a length (they do not collapse to reach the angles)
    assert min(formdiagram.edge_length(edge) for edge in formdiagram.edges()) > 0.5
    # weight 1 changes the form diagram only
    assert all(forcediagram.vertex_coordinates(vertex) == xyz for vertex, xyz in force_xyz.items())


def test_volmesh_reciprocate_weight_0_perturbed_force(forcediagram, formdiagram):
    # only the interior vertices move: the boundary faces keep their initial normals during
    # reciprocation, so they stay planar and the meshgrid is a solution
    rng = random.Random(1)
    for vertex in sorted(forcediagram.vertices()):
        if forcediagram.is_vertex_on_boundary(vertex):
            continue
        x, y, z = forcediagram.vertex_coordinates(vertex)
        forcediagram.vertex_attributes(vertex, "xyz", [x + rng.uniform(-0.3, 0.3), y + rng.uniform(-0.3, 0.3), z + rng.uniform(-0.3, 0.3)])
    form_xyz = {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()}
    boundary_normals = {halfface: forcediagram.halfface_normal(halfface) for halfface in forcediagram.halffaces_on_boundaries()}

    before = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )
    volmesh_reciprocate(forcediagram, formdiagram, kmax=100, weight=0.0, tolerance=TOLERANCE, tolerance_boundary=TOLERANCE_BOUNDARY)
    after = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )

    assert before > 10
    assert after < TOLERANCE
    assert max(angle_vectors(forcediagram.halfface_normal(halfface), normal, deg=True) for halfface, normal in boundary_normals.items()) < TOLERANCE_BOUNDARY
    for cell in forcediagram.cells():
        mesh = forcediagram.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0
    # weight 0 changes the force diagram only
    assert all(formdiagram.node_coordinates(node) == xyz for node, xyz in form_xyz.items())


@pytest.mark.parametrize("perturbed", ["form", "force"])
def test_volmesh_reciprocate_weight_half_lowers_deviation(forcediagram, formdiagram, perturbed):
    rng = random.Random(1)
    if perturbed == "form":
        for node in sorted(formdiagram.nodes()):
            x, y, z = formdiagram.node_coordinates(node)
            formdiagram.node_attributes(node, "xyz", [x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), z + rng.uniform(-0.5, 0.5)])
    else:
        for vertex in sorted(forcediagram.vertices()):
            x, y, z = forcediagram.vertex_coordinates(vertex)
            forcediagram.vertex_attributes(vertex, "xyz", [x + rng.uniform(-0.3, 0.3), y + rng.uniform(-0.3, 0.3), z + rng.uniform(-0.3, 0.3)])
    force_xyz = {vertex: forcediagram.vertex_coordinates(vertex) for vertex in forcediagram.vertices()}
    form_xyz = {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()}

    before = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )
    volmesh_reciprocate(forcediagram, formdiagram, kmax=10, weight=0.5, tolerance=TOLERANCE)
    after = max(
        angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges()
    )

    assert after < before
    # weight 0.5 changes both diagrams
    assert max(distance_point_point(forcediagram.vertex_coordinates(vertex), xyz) for vertex, xyz in force_xyz.items()) > TOL.absolute
    assert max(distance_point_point(formdiagram.node_coordinates(node), xyz) for node, xyz in form_xyz.items()) > TOL.absolute


def test_volmesh_reciprocate_fix_vkeys(forcediagram, formdiagram):
    rng = random.Random(1)
    for node in sorted(formdiagram.nodes()):
        x, y, z = formdiagram.node_coordinates(node)
        formdiagram.node_attributes(node, "xyz", [x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), z + rng.uniform(-0.5, 0.5)])
    # the node of the middle cell
    fixed = formdiagram.node_coordinates(13)

    volmesh_reciprocate(forcediagram, formdiagram, kmax=100, weight=1.0, fix_vkeys=[13], tolerance=TOLERANCE)

    assert (
        max(angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True) for u, v in formdiagram.edges())
        < TOLERANCE
    )
    assert formdiagram.node_coordinates(13) == fixed


def test_volmesh_reciprocate_callback(forcediagram, formdiagram):
    rng = random.Random(1)
    for node in sorted(formdiagram.nodes()):
        x, y, z = formdiagram.node_coordinates(node)
        formdiagram.node_attributes(node, "xyz", [x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), z + rng.uniform(-0.5, 0.5)])
    calls = []
    args = {"name": "args"}

    volmesh_reciprocate(
        forcediagram,
        formdiagram,
        kmax=5,
        weight=1.0,
        tolerance=TOLERANCE,
        callback=lambda volmesh, network, k, callback_args: calls.append((volmesh, network, k, callback_args)),
        callback_args=args,
    )

    assert [k for _, _, k, _ in calls] == [0, 1, 2, 3, 4]
    assert all(volmesh is forcediagram and network is formdiagram and callback_args is args for volmesh, network, _, callback_args in calls)


def test_volmesh_reciprocate_callback_not_callable_raises(forcediagram, formdiagram):
    with pytest.raises(Exception, match="not callable"):
        volmesh_reciprocate(forcediagram, formdiagram, kmax=1, callback=1)
