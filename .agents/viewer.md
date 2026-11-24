# Viewer: drawing models in compas_viewer

`from compas_3gs.viewer import GSViewer` — a `compas_viewer.Viewer`. Only this module needs
`compas_viewer`; the rest of the package does not.

```python
viewer = GSViewer()
groups = viewer.add_model(model)  # form diagram + force cells (form first) or force diagram (force first)
viewer.show()
```

| Call | |
|---|---|
| `add_model(model, labels=True, offset=None)` | everything; `offset` places the force diagram (default: beside the form along x) |
| `add_form(form, labels=True, offset=None)` | groups `form/members, labels, supports, loads, reactions`; members red in compression, blue in tension, width by force |
| `update_form(form, groups)` | after the forces change (for example in a slider callback): recolours and resizes the members and refreshes their info; labels, loads and reactions keep their drawn state |
| `add_cells(model, offset)` | the force polyhedron of every node of `model.cells`, scaled together, faces coloured like their members (load black, reaction green) |
| `add_force_diagram(force, form, offset)` | a `ForceVolMesh`, one mesh per cell |
| `add_support(point, dofs, size, parent=None)` | support glyph: tetrahedron; + two rollers for one restrained direction, + a line for two |
| `member_info(form, edge)` | `{"type", "force", "fx", "fy", "fz"}` for the info panel |
| `viewer.info`, `show_info(sceneform, node)` | the info panel of the sidebar: shows `obj.attributes["data"]` of the selected object (tree or view) |

Every `add_*` returns its groups by name, so parts can be hidden (`groups["reactions"].show = False`).
Any scene object can fill the info panel by setting `obj.attributes["data"] = {...}` (one flat dict).

A live slider (see `examples/model/05_warren_truss_fem.py`):

```python
from compas_viewer.components import Slider

def pick(slider, value):
    model.solve(reactions={(node, "z"): value})
    viewer.update_form(model.form, groups)

viewer.ui.sidedock.show = True
viewer.ui.sidedock.add(Slider(title="Bearing (kN)", value=0.0, min_val=0.0, max_val=60.0, step=0.1, action=pick))
```

## compas_viewer facts (2.0.2)

- `Viewer` is a singleton keyed by `Viewer`: create one viewer per script. Creating `Viewer()` and
  then `GSViewer()` in one process returns the first one (tests reset `SingletonMeta._instances`).
- No VolMesh scene object: draw cells with `volmesh.cell_to_mesh(cell)`.
- `GraphObject` ignores per-edge colours: draw members as `Line` objects (as `add_form` does).
- Colours: `Color` objects, or 0–255 tuples; `show_lines=`, not `show_edges=`.
- Animation: `@viewer.on(interval=..., frames=...)` and `obj.update(update_data=True)`.
- Tag text and vector anchors cannot change after drawing.
- Headless runs: `QT_QPA_PLATFORM=offscreen` and replace `Viewer.show` (see `tests/test_examples.py`).
