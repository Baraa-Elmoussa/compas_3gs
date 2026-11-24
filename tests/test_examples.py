import pathlib
import runpy

import pytest

compas_viewer = pytest.importorskip("compas_viewer", reason="compas_viewer is an optional dependency, needed by the examples")

ROOT = pathlib.Path(__file__).resolve().parent.parent
FOLDERS = ["cell_manipulation", "multi_cell_polyhedron", "physical_structures"]
EXAMPLES = sorted(path for folder in FOLDERS for path in (ROOT / "examples" / folder).glob("*.py"))


@pytest.mark.gui
@pytest.mark.parametrize("path", EXAMPLES, ids=[path.relative_to(ROOT / "examples").as_posix() for path in EXAMPLES])
def test_example(path, monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    monkeypatch.setattr(compas_viewer.Viewer, "show", lambda self: None)
    # the examples are run from the repo root, as their docstrings say
    monkeypatch.chdir(ROOT)

    runpy.run_path(str(path), run_name="__main__")
