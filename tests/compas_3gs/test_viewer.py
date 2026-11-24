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
        assert member.linecolor == Color.red()  # all legs are in compression


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
    viewer.add_support([0.0, 0.0, 0.0], dofs, 0.1, parent=group)
    assert len(group.children) == objects
