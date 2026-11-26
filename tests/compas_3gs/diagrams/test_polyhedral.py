import random
from collections import Counter

import pytest

import compas
from compas.geometry import angle_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.datastructures import Network3gs
from compas_3gs.datastructures import VolMesh3gs
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork


@pytest.fixture
def force():
    # 2 x 2 x 2 cubic cells of size 0.5
    return ForceVolMesh.from_meshgrid(1, 1, 1, 2, 2, 2)


@pytest.fixture
def form(force):
    # one node per cell at its centroid, one edge per pair of cells that share a face
    form = FormNetwork()
    for cell in force.cells():
        x, y, z = force.cell_centroid(cell)
        form.add_node(cell, x=x, y=y, z=z)
    for cell in force.cells():
        for nbr in dict.fromkeys(force.cell_neighbors(cell)):
            if cell < nbr:
                form.add_edge(cell, nbr)
    form.dual = force
    force.primal = form
    return form


# ==============================================================================
# FormNetwork
# ==============================================================================


def test_formnetwork_constructor():
    form = FormNetwork()
    assert isinstance(form, Network3gs)
    assert form.name == "FormNetwork"
    assert form.dual is None
    assert form.attributes == {"scale": 1}
    assert form.default_node_attributes == {
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
        "is_fixed": False,
        "x_fix": False,
        "y_fix": False,
        "z_fix": False,
        "support": "",
        "px": 0.0,
        "py": 0.0,
        "pz": 0.0,
        "_rx": 0.0,
        "_ry": 0.0,
        "_rz": 0.0,
    }
    for name, value in {"l_min": 1, "l_max": 20, "target_vector": None, "target_length": None, "_a": 0.0, "_f": 0.0}.items():
        assert form.default_edge_attributes[name] == value
    assert FormNetwork(name="form").name == "form"


def test_formnetwork_vertex_update_xyz(form):
    form.vertex_update_xyz(0, [1.0, 2.0, 3.0])
    assert TOL.is_allclose(form.node_coordinates(0), [1.0, 2.0, 3.0])
    form.node_attribute(1, "y_fix", True)
    form.vertex_update_xyz(1, [1.0, 2.0, 3.0])
    assert TOL.is_allclose(form.node_coordinates(1), [1.0, 0.75, 3.0])


def test_dual_halfface(form, force):
    assert form.number_of_edges() == 12
    for u, v in form.edges():
        halfface = form.dual_halfface((u, v))
        assert halfface in force.cell_faces(u)
        assert force.halfface_opposite_cell(halfface) == v
        assert halfface == force.cell_pair_halffaces(u, v)[0]
        assert form.dual_halfface((v, u)) == force.halfface_opposite_halfface(halfface)
    with pytest.raises(KeyError):
        form.dual_halfface((0, 7))


def test_update_angle_deviations_on_a_reciprocal_pair(form, force):
    form.update_angle_deviations()
    for edge in form.edges():
        halfface = form.dual_halfface(edge)
        assert TOL.is_zero(form.edge_attribute(edge, "_a"))
        assert TOL.is_zero(force.face_attribute(halfface, "_a"))
        assert TOL.is_zero(force.face_attribute(force.halfface_opposite_halfface(halfface), "_a"))


def test_update_angle_deviations_on_a_perturbed_form(form, force):
    rng = random.Random(1)
    for node in form.nodes():
        form.node_attributes(node, "xyz", [axis + rng.uniform(-0.05, 0.05) for axis in form.node_coordinates(node)])
    form.update_angle_deviations()
    for edge in form.edges():
        halfface = form.dual_halfface(edge)
        angle = angle_vectors(form.edge_direction(edge), force.face_normal(halfface), deg=True)
        assert angle > 1.0
        assert form.edge_attribute(edge, "_a") == pytest.approx(angle)
        assert force.face_attribute(halfface, "_a") == pytest.approx(angle)


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_formnetwork_roundtrip(form, how):
    form.name = "grid form"
    form.attributes["scale"] = 2.0
    form.node_attribute(0, "is_fixed", True)
    form.edge_attribute((0, 1), "target_length", 0.75)
    form.edge_attribute((0, 2), "target_vector", [1.0, 0.0, 0.0])
    form.set_support(0, "slider")
    form.add_force(7, [0.0, 0.0, -1.0])
    form.edge_attribute((0, 1), "_f", -2.5)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(form))
    else:
        other = form.copy(copy_guid=True)

    assert type(other) is FormNetwork
    assert other.guid == form.guid
    assert other.name == "grid form"
    assert other.dual is None
    assert other.attributes == {"scale": 2.0}
    assert other.default_node_attributes == form.default_node_attributes
    assert other.default_edge_attributes == form.default_edge_attributes
    assert list(other.nodes()) == list(form.nodes())
    assert list(other.edges()) == list(form.edges())
    for node in form.nodes():
        assert TOL.is_allclose(other.node_coordinates(node), form.node_coordinates(node))
        assert other.node_attributes(node) == form.node_attributes(node)
    for edge in form.edges():
        assert other.edge_attributes(edge) == form.edge_attributes(edge)


def test_node_at(form):
    point = form.node_coordinates(5)
    assert form.node_at(point) == 5
    with pytest.raises(KeyError):
        form.node_at([10.0, 10.0, 10.0])
    with pytest.raises(KeyError):
        form.node_at([point[0] + 0.01, point[1], point[2]])
    assert form.node_at([point[0] + 0.01, point[1], point[2]], tol=0.1) == 5


@pytest.mark.parametrize(
    "dofs, stored",
    [("pinned", "xyz"), ("slider", "xz"), ("roller", "z"), ("free", ""), ("zx", "xz"), ("y", "y"), ("xyz", "xyz")],
)
def test_set_support(form, dofs, stored):
    form.set_support(3, dofs)
    assert form.node_attribute(3, "support") == stored
    assert list(form.supports()) == ([3] if stored else [])


def test_set_support_rejects_other_names(form):
    with pytest.raises(ValueError):
        form.set_support(3, "fixed")


def test_add_force_adds_up(form):
    form.add_force(2, [1.0, 0.0, -2.0])
    form.add_force(2, [0.5, 1.0, 0.0])
    assert TOL.is_allclose(form.node_load(2), [1.5, 1.0, -2.0])
    assert TOL.is_allclose(form.node_load(3), [0.0, 0.0, 0.0])


def test_node_forces_directions_and_signs(form):
    # node 0 is the corner cell of the grid; its neighbours are 1, 2 and 4
    form.edge_attribute((0, 1), "_f", 2.0)  # tension pulls node 0 towards node 1
    form.edge_attribute((0, 2), "_f", -3.0)  # compression pushes node 0 away from node 2
    form.add_force(0, [0.0, 0.0, -1.0])
    form.node_attributes(0, ["_rx", "_ry", "_rz"], [0.0, 0.0, 1.0])

    forces = form.node_forces(0)

    assert set(forces) == {(0, 1), (0, 2), (0, 4), "load", "reaction"}
    towards_1 = [b - a for a, b in zip(form.node_coordinates(0), form.node_coordinates(1))]
    towards_2 = [b - a for a, b in zip(form.node_coordinates(0), form.node_coordinates(2))]
    assert TOL.is_allclose(forces[(0, 1)], [2.0 * c / 0.5 for c in towards_1])
    assert TOL.is_allclose(forces[(0, 2)], [-3.0 * c / 0.5 for c in towards_2])
    assert TOL.is_zero(forces[(0, 4)].length)
    assert TOL.is_allclose(forces["load"], [0.0, 0.0, -1.0])
    assert TOL.is_allclose(forces["reaction"], [0.0, 0.0, 1.0])
    assert "load" not in form.node_forces(3)


# ==============================================================================
# ForceVolMesh
# ==============================================================================


def test_forcevolmesh_constructor():
    force = ForceVolMesh()
    assert isinstance(force, VolMesh3gs)
    assert force.name == "ForceVolMesh"
    assert force.primal is None
    assert force.attributes == {"scale": 1, "convention": -1}
    assert force.default_vertex_attributes == {"x": 0.0, "y": 0.0, "z": 0.0, "x_fix": False, "y_fix": False, "z_fix": False, "is_fixed": False}
    assert force.default_edge_attributes == {"l_min": 1, "l_max": 20, "target_vector": None, "target_length": None}
    assert force.default_face_attributes == {"_a": 0.0, "target_area": None, "target_normal": None}
    assert force.default_cell_attributes == {"dir": None}
    assert ForceVolMesh(name="force").name == "force"


def test_forcevolmesh_scale_and_convention_are_attributes(force):
    assert force.attributes["scale"] == 1
    assert force.attributes["convention"] == -1
    force.attributes["scale"] = 0.5
    force.attributes["convention"] = 1
    assert force.attributes == {"scale": 0.5, "convention": 1}
    assert ForceVolMesh().attributes == {"scale": 1, "convention": -1}


def test_forcevolmesh_from_vertices_and_cells_is_valid(force):
    rng = random.Random(2)
    vertices = [[axis + rng.uniform(-0.1, 0.1) for axis in force.vertex_coordinates(vertex)] for vertex in force.vertices()]
    cells = [[force.halfface_vertices(halfface) for halfface in force.cell_faces(cell)] for cell in force.cells()]
    volmesh = ForceVolMesh.from_vertices_and_cells(vertices, cells)

    assert type(volmesh) is ForceVolMesh
    assert volmesh.attributes == {"scale": 1, "convention": -1}
    assert volmesh.cell_attribute(0, "dir") is None
    assert volmesh.face_attribute(0, "_a") == 0.0
    assert volmesh.number_of_cells() == 8
    assert len(list(volmesh.halffaces())) == 48
    assert len(volmesh.halffaces_on_boundaries()) == 24
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0
    shared = Counter(frozenset(volmesh.halfface_vertices(halfface)) for halfface in volmesh.halffaces())
    for halfface in volmesh.halffaces():
        opposite = volmesh.halfface_opposite_halfface(halfface)
        if shared[frozenset(volmesh.halfface_vertices(halfface))] == 1:
            assert opposite is None
        else:
            assert opposite is not None
            assert volmesh.halfface_opposite_halfface(opposite) == halfface


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_forcevolmesh_roundtrip(force, how):
    force.name = "grid force"
    force.attributes["scale"] = 2.0
    force.vertex_attribute(0, "is_fixed", True)
    force.face_attribute(0, "target_area", 0.5)
    force.cell_attribute(0, "dir", 1)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(force))
    else:
        other = force.copy(copy_guid=True)

    assert type(other) is ForceVolMesh
    assert other.guid == force.guid
    assert other.name == "grid force"
    assert other.primal is None
    assert other.attributes == {"scale": 2.0, "convention": -1}
    assert other.default_vertex_attributes == force.default_vertex_attributes
    assert other.default_edge_attributes == force.default_edge_attributes
    assert other.default_face_attributes == force.default_face_attributes
    assert other.default_cell_attributes == force.default_cell_attributes
    assert list(other.vertices()) == list(force.vertices())
    assert list(other.cells()) == list(force.cells())
    assert list(other.halffaces()) == list(force.halffaces())
    for vertex in force.vertices():
        assert TOL.is_allclose(other.vertex_coordinates(vertex), force.vertex_coordinates(vertex))
        assert other.vertex_attributes(vertex) == force.vertex_attributes(vertex)
    for cell in force.cells():
        assert other.cell_faces(cell) == force.cell_faces(cell)
        assert other.cell_attributes(cell) == force.cell_attributes(cell)
    for halfface in force.halffaces():
        assert other.halfface_vertices(halfface) == force.halfface_vertices(halfface)
        assert other.face_attributes(halfface) == force.face_attributes(halfface)
