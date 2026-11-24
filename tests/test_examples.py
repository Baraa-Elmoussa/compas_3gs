import pathlib
import runpy

import pytest

compas_viewer = pytest.importorskip("compas_viewer", reason="compas_viewer is an optional dependency, needed by the examples")

from compas_viewer.singleton import SingletonMeta  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
FOLDERS = ["model", "cell_manipulation", "multi_cell_polyhedron", "physical_structures"]
EXAMPLES = sorted(path for folder in FOLDERS for path in (ROOT / "examples" / folder).glob("*.py"))


@pytest.mark.gui
@pytest.mark.parametrize("path", EXAMPLES, ids=[path.relative_to(ROOT / "examples").as_posix() for path in EXAMPLES])
def test_example(path, monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    monkeypatch.setattr(compas_viewer.Viewer, "show", lambda self: None)
    # compas_viewer keeps one viewer per process, keyed by the first subclass of its Singleton (Viewer), and one
    # QApplication: give every example a new viewer of its own class (Viewer or GSViewer) on the same application
    monkeypatch.setattr(SingletonMeta, "_instances", {})
    app = QApplication.instance()
    if app is not None:
        monkeypatch.setattr(compas_viewer.Viewer, "create_app", lambda self: app)
    # the examples are run from the repo root, as their docstrings say
    monkeypatch.chdir(ROOT)

    runpy.run_path(str(path), run_name="__main__")
