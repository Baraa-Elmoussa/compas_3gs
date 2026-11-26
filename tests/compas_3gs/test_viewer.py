from math import cos
from math import radians
from math import sin

import pytest

compas_viewer = pytest.importorskip("compas_viewer", reason="compas_viewer is an optional dependency, needed by compas_3gs.viewer")

from compas_viewer.scene.lineobject import LineObject  # noqa: E402
from compas_viewer.singleton import SingletonMeta  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from compas.colors import Color  # noqa: E402
from compas_3gs.diagrams import ForceVolMesh  # noqa: E402
from compas_3gs.model import GSModel  # noqa: E402
from compas_3gs.operations import volmesh_cell_subdivide_barycentric  # noqa: E402
from compas_3gs.viewer import GSViewer  # noqa: E402

APEX = [0.0, 0.0, 3.0]
FEET = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (0, 120, 240)]


@pytest.fixture
def viewer(monkeypatch):
    # compas_viewer keeps one viewer per process, keyed by the first subclass of its Singleton (Viewer), and one
    # QApplication: give every test a new GSViewer on the same application
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    monkeypatch.setattr(SingletonMeta, "_instances", {})
    app = QApplication.instance()
    if app is not None:
        monkeypatch.setattr(compas_viewer.Viewer, "create_app", lambda self: app)
    return GSViewer()


@pytest.fixture
def tripod():
    model = GSModel.from_lines([(APEX, foot) for foot in FEET])
    for foot in FEET:
        model.set_support(model.node_at(foot), "pinned")
    model.add_force(model.node_at(APEX), [1.5, 0.0, -10.0])
    model.solve()
    model.derive_force_diagram()
    return model


@pytest.mark.gui
def test_add_model_form_first(viewer, tripod):
    groups = viewer.add_model(tripod)

    assert isinstance(viewer, GSViewer)
    assert set(groups) == {"form", "members", "labels", "supports", "loads", "reactions", "force", "cells"}
    assert len(groups["members"].children) == 3
    assert len(groups["labels"].children) == 3
    assert len(groups["supports"].children) == 3  # one tetrahedron per pinned support
    assert len(groups["loads"].children) == 1
    assert len(groups["reactions"].children) == 3
    assert len(groups["cells"].children) == 1  # the feet have no force polyhedron
    for member in groups["members"].children:
        assert isinstance(member, LineObject)
        assert member.linecolor == Color.blue()  # all legs are in compression


@pytest.mark.gui
def test_add_model_without_labels(viewer, tripod):
    groups = viewer.add_model(tripod, labels=False)
    assert len(groups["labels"].children) == 0


@pytest.mark.gui
def test_add_model_force_first(viewer):
    force = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
    volmesh_cell_subdivide_barycentric(force, next(iter(force.cells())))
    model = GSModel.from_force_diagram(force)

    groups = viewer.add_model(model)

    assert len(groups["members"].children) == 12
    assert len(groups["loads"].children) == 6
    assert len(groups["supports"].children) == 0
    assert len(groups["cells"].children) == 6  # one mesh per cell of the force diagram


@pytest.mark.gui
@pytest.mark.parametrize("dofs, objects", [("xyz", 1), ("xz", 2), ("z", 3)], ids=["pinned", "slider", "roller"])
def test_support_glyphs(viewer, dofs, objects):
    # a tetrahedron, with a line along the free direction for a slider and two rollers for a roller
    group = viewer.scene.add_group(name="supports")
    glyph = viewer.add_support([0.0, 0.0, 0.0], dofs, 0.1, parent=group)
    assert len(group.children) == objects
    assert glyph == list(group.children)


@pytest.mark.gui
def test_members_carry_their_info(viewer, tripod):
    groups = viewer.add_model(tripod)

    for member in groups["members"].children:
        edge = member.attributes["edge"]
        data = member.attributes["data"]
        assert list(data) == ["type", "force", "fx", "fy", "fz"]
        force = tripod.form.edge_force(edge)
        assert data["type"] == "compression"
        assert data["force"] == pytest.approx(force, abs=1e-3)
        # the force the member applies to its first node, along the edge direction
        direction = tripod.form.edge_direction(edge)
        assert [data["fx"], data["fy"], data["fz"]] == pytest.approx([force * c for c in direction], abs=1e-3)


@pytest.mark.gui
def test_node_info(viewer, tripod):
    form = tripod.form
    apex = tripod.node_at(APEX)

    info = viewer.node_info(form, apex)

    assert list(info)[:4] == ["node", "px", "py", "pz"]
    assert [info["px"], info["py"], info["pz"]] == pytest.approx([1.5, 0.0, -10.0])
    assert "support" not in info
    assert {key: value for key, value in info.items() if key.startswith("member")} == pytest.approx(
        {"member {}-{}".format(u, v): form.edge_force((u, v)) for u, v in form.node_edges(apex)}, abs=1e-3
    )
    # the reactions of the three feet balance the load
    feet = [viewer.node_info(form, tripod.node_at(foot)) for foot in FEET]
    for foot in feet:
        assert list(foot)[:5] == ["node", "support", "rx", "ry", "rz"]
        assert foot["support"] == form.node_attribute(foot["node"], "support")
    assert [sum(foot[key] for foot in feet) for key in ("rx", "ry", "rz")] == pytest.approx([-1.5, 0.0, 10.0], abs=1e-2)


@pytest.mark.gui
def test_loads_supports_reactions_and_cells_carry_the_forces_on_their_node(viewer, tripod):
    groups = viewer.add_model(tripod)

    for name in ("loads", "supports", "reactions", "cells"):
        for obj in groups[name].children:
            assert obj.attributes["data"] == viewer.node_info(tripod.form, obj.attributes["node"])
    assert [obj.name for obj in groups["loads"].children] == ["load {}".format(tripod.node_at(APEX))]
    assert sorted(obj.name for obj in groups["reactions"].children) == sorted("reaction {}".format(tripod.node_at(foot)) for foot in FEET)
    assert groups["cells"].children[0].attributes["node"] == tripod.node_at(APEX)


@pytest.mark.gui
def test_force_diagram_cells_carry_the_forces_on_their_node(viewer):
    force = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
    volmesh_cell_subdivide_barycentric(force, next(iter(force.cells())))
    model = GSModel.from_force_diagram(force)

    groups = viewer.add_model(model)

    for obj in groups["cells"].children:
        node = obj.attributes["node"]
        assert obj.attributes["data"] == viewer.node_info(model.form, node)
        assert sum(key.startswith("member") for key in obj.attributes["data"]) == model.form.degree(node)


@pytest.mark.gui
def test_info_panel_shows_the_selected_member(viewer, tripod):
    groups = viewer.add_model(tripod)
    member = groups["members"].children[0]

    viewer.show_info(viewer.ui.sidebar.sceneform, member)

    shown = [(node.name, node.attributes.get("value")) for node in viewer.info.tree.nodes if not node.is_root]
    assert shown == list(member.attributes["data"].items())
    assert viewer.ui.sidebar.sceneform.action == viewer.show_info


@pytest.mark.gui
def test_update_form_after_solving_again(viewer, tripod):
    groups = viewer.add_form(tripod.form, labels=False)
    tripod.form.node_attributes(tripod.node_at(APEX), ["px", "py", "pz"], [0.0, 0.0, 10.0])  # pull the apex up: tension
    tripod.solve()

    viewer.update_form(tripod.form, groups)

    for member in groups["members"].children:
        assert member.attributes["data"]["type"] == "tension"
        assert member.linecolor == Color.red()
    for name in ("loads", "supports", "reactions"):
        for obj in groups[name].children:
            assert obj.attributes["data"] == viewer.node_info(tripod.form, obj.attributes["node"])
    assert sum(obj.attributes["data"]["rz"] for obj in groups["reactions"].children) == pytest.approx(-10.0, abs=1e-2)


@pytest.mark.gui
def test_add_form_with_an_offset(viewer, tripod):
    groups = viewer.add_form(tripod.form, offset=[0.0, -5.0, 1.0])
    member = groups["members"].children[0]
    u, v = member.attributes["edge"]
    assert list(member.geometry.start) == pytest.approx([a + b for a, b in zip(tripod.form.node_coordinates(u), [0.0, -5.0, 1.0])])
