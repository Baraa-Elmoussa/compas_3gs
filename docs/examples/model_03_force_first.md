# Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces

A cube split into six pyramids is the force diagram of an octahedral frame under six loads.
`GSModel.from_force_diagram` makes the dual network the form diagram, reciprocates it so that every
member is perpendicular to its face, takes the member forces from the face areas and the loads
from the boundary faces.

![Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces](../assets/images/examples/model_03_force_first.png)

```python title="examples/model/03_force_first.py"
--8<-- "examples/model/03_force_first.py"
```
