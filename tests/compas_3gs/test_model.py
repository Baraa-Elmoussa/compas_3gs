from math import cos
from math import radians
from math import sin

import pytest

import compas
from compas.geometry import angle_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.algorithms import EquilibriumResult
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.model import GSModel
from compas_3gs.operations import volmesh_cell_subdivide_barycentric

APEX = [0.0, 0.0, 3.0]
FEET = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (0, 120, 240)]
LOAD = [1.5, 0.0, -10.0]


@pytest.fixture
def tripod():
    model = GSModel.from_lines([(APEX, foot) for foot in FEET])
    for foot in FEET:
        model.set_support(node=model.node_at(foot), dofs="pinned")
    model.add_force(node=model.node_at(APEX), vector=LOAD)
    return model


@pytest.fixture
def force_first():
    force = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
    volmesh_cell_subdivide_barycentric(force, next(iter(force.cells())))
    return GSModel.from_force_diagram(force)


def test_constructor():
    model = GSModel()
    assert isinstance(model.form, FormNetwork)
    assert model.force is None
    assert model.result is None
    assert model.cells == {} and model.cell_faces == {} and model.skipped == {}


def test_solve_stores_forces_and_reactions(tripod):
    result = tripod.solve()

    assert isinstance(result, EquilibriumResult)
    assert tripod.result is result
    assert result.determinate
    forces = tripod.forces()
    assert len(forces) == 3
    assert all(force < 0 for force in forces.values())
    for edge, force in forces.items():
        assert tripod.form.edge_attribute(edge, "_f") == force
    reactions = tripod.reactions()
    assert set(reactions) == {tripod.node_at(foot) for foot in FEET}
    total = [sum(reaction[k] for reaction in reactions.values()) + LOAD[k] for k in range(3)]
    assert TOL.is_allclose(total, [0, 0, 0])


def test_every_node_is_in_equilibrium_after_solve(tripod):
    tripod.solve()
    for node in tripod.form.nodes():
        total = [sum(force[k] for force in tripod.form.node_forces(node).values()) for k in range(3)]
        assert TOL.is_allclose(total, [0, 0, 0])


def test_derive_force_diagram(tripod):
    tripod.solve()
    cells = tripod.derive_force_diagram()

    apex = tripod.node_at(APEX)
    assert set(cells) == {apex}
    assert tripod.cells is cells
    assert set(tripod.skipped) == {tripod.node_at(foot) for foot in FEET}
    assert all("four or more" in reason for reason in tripod.skipped.values())
    cell = cells[apex]
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    assert sorted(map(str, tripod.cell_faces[apex].values())) == sorted(map(str, list(tripod.form.edges()) + ["load"]))
    forces = tripod.form.node_forces(apex)
    for face, name in tripod.cell_faces[apex].items():
        assert cell.face_area(face) == pytest.approx(forces[name].length, rel=1e-3)
        assert angle_vectors(cell.face_normal(face), forces[name], deg=True) < 0.1
    check = tripod.check()
    assert set(check) == {apex}
    assert check[apex]["area_error"] < 1e-3
    assert check[apex]["angle"] < 0.1


def test_solve_removes_old_cells(tripod):
    tripod.solve()
    tripod.derive_force_diagram()
    tripod.solve()
    assert tripod.cells == {} and tripod.cell_faces == {} and tripod.skipped == {}


def test_supports_of_every_type():
    model = GSModel.from_lines([([0, 0, 0], [1, 0, 0]), ([1, 0, 0], [0, 1, 0]), ([0, 1, 0], [0, 0, 0]), ([0, 0, 1], [0, 0, 0]), ([0, 0, 1], [1, 0, 0]), ([0, 0, 1], [0, 1, 0])])
    model.set_support(model.node_at([0, 0, 0]), "pinned")
    model.set_support(model.node_at([1, 0, 0]), "yz")  # a slider along the base member
    model.set_support(model.node_at([0, 1, 0]), "roller")
    model.add_force(model.node_at([0, 0, 1]), [0.3, 0.2, -1.0])

    result = model.solve()

    assert result.determinate
    reactions = model.reactions()
    assert reactions[model.node_at([1, 0, 0])][0] == 0.0
    assert reactions[model.node_at([0, 1, 0])][0] == 0.0 and reactions[model.node_at([0, 1, 0])][1] == 0.0
    total = [sum(reaction[k] for reaction in reactions.values()) for k in range(3)]
    assert TOL.is_allclose(total, [-0.3, -0.2, 1.0])


def test_solve_raises_for_an_indeterminate_truss():
    corners = [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0]]
    model = GSModel.from_lines([([0.0, 0.0, 2.0], corner) for corner in corners])
    for corner in corners:
        model.set_support(model.node_at(corner), "pinned")
    model.add_force(model.node_at([0.0, 0.0, 2.0]), [0.0, 0.0, -10.0])
    with pytest.raises(ValueError, match="indeterminate"):
        model.solve()


def test_solve_raises_for_a_mechanism():
    model = GSModel.from_lines([([0.0, 0.0, 2.0], [1.0, 1.0, 0.0]), ([0.0, 0.0, 2.0], [-1.0, -1.0, 0.0])])
    model.set_support(model.node_at([1.0, 1.0, 0.0]), "pinned")
    model.set_support(model.node_at([-1.0, -1.0, 0.0]), "pinned")
    model.add_force(model.node_at([0.0, 0.0, 2.0]), [1.0, -1.0, -10.0])
    with pytest.raises(ValueError, match="mechanism"):
        model.solve()


def test_from_force_diagram(force_first):
    form = force_first.form
    assert isinstance(form, FormNetwork)
    assert form.dual is force_first.force
    assert force_first.force.primal is form
    assert form.number_of_nodes() == 6
    assert form.number_of_edges() == 12
    # every member is perpendicular to its face and in compression, and every node is in equilibrium
    form.update_angle_deviations()
    assert max(form.edge_attribute(edge, "_a") for edge in form.edges()) < 0.01
    assert all(force < 0 for force in force_first.forces().values())
    for node in form.nodes():
        total = [sum(force[k] for force in form.node_forces(node).values()) for k in range(3)]
        assert TOL.is_allclose(total, [0, 0, 0])


def test_summary(tripod):
    tripod.solve()
    summary = tripod.summary()
    assert summary.count(" C") == 3
    assert "reactions:" in summary


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(tripod, how):
    tripod.name = "tripod"
    tripod.solve()

    other = compas.json_loads(compas.json_dumps(tripod)) if how == "json" else tripod.copy(copy_guid=True)

    assert type(other) is GSModel
    assert other.guid == tripod.guid
    assert other.name == "tripod"
    assert other.force is None
    assert list(other.form.nodes()) == list(tripod.form.nodes())
    assert list(other.form.edges()) == list(tripod.form.edges())
    for node in tripod.form.nodes():
        assert other.form.node_attributes(node) == tripod.form.node_attributes(node)
    assert other.forces() == tripod.forces()


@pytest.mark.roundtrip
@pytest.mark.xfail(strict=True, reason="compas 2.15.1 VolMesh.__data__ writes each face of a cell once per halfedge, so the force diagram comes back with extra faces")
def test_roundtrip_force_first(force_first):
    other = compas.json_loads(compas.json_dumps(force_first))
    assert other.form.dual is other.force
    assert len(list(other.force.halffaces())) == len(list(force_first.force.halffaces()))
