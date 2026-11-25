# Viewing

`GSViewer` is a compas_viewer `Viewer` that draws a model. It is the one part of the package that needs
compas_viewer (`pip install compas_viewer`).
{ .lead }

```python
from compas_3gs.viewer import GSViewer

viewer = GSViewer()
groups = viewer.add_model(model)
viewer.show()
```

## What is drawn

Every part goes into a group of the scene tree, so it can be shown or hidden on its own:

- `form`
    - `members`: red in compression, blue in tension, grey without force, thicker with a larger force;
    - `labels`: the member forces (`add_model(model, labels=False)` leaves them out);
    - `supports`: a tetrahedron under each support, on two rollers if it restrains one direction and on a
      line along the free direction if it restrains two;
    - `loads`: black arrows; `reactions`: green arrows.
- `force`
    - `cells`: the force polyhedra of a form-first model, beside the form diagram, faces coloured like their
      members, or the polyhedral force diagram of a force-first model.

`add_model` returns the groups by name. `add_form(form, offset=...)`, `add_cells(model, offset)`,
`add_force_diagram(force, form, offset)` and `add_support(point, dofs, size)` draw the parts one by one.

## The info panel

Selecting a member, in the scene tree or in the view, shows in the sidebar:

| Name | Value |
|---|---|
| `type` | `compression`, `tension` or `zero` |
| `force` | the axial force, tension positive |
| `fx`, `fy`, `fz` | the force the member applies to its first node |

The panel shows the attribute `data` of any selected scene object, so a script can fill it for its own
objects: `obj.attributes["data"] = {"name": value, ...}`.

## Changing the forces while viewing

After solving again, `viewer.update_form(form, groups)` recolours and resizes the members and refreshes
their info. With a slider of compas_viewer, a script can pick the state of an indeterminate structure
while it is shown:

```python
from compas_viewer.components import Slider

def pick(slider, value):
    model.solve(reactions={(node, "z"): value})
    viewer.update_form(model.form, groups)

viewer.ui.sidedock.show = True
viewer.ui.sidedock.add(Slider(title="Bearing (kN)", value=0.0, min_val=0.0, max_val=60.0, step=0.1, action=pick))
```

Labels and arrows keep the state in which they were drawn. The full script is
[Warren truss, FEM and 3GS](../examples/model_05_warren_truss_fem.md).
