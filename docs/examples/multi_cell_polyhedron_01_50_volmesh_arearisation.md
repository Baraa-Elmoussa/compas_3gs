# Arearise chosen halffaces of a volmesh to a target area

The top and bottom halffaces of the centre cell of a 3x3 mat are resized to a target area while
all face normals are kept. The outer boundary does not move, so the other cells change size to make
room (the corner cells shrink, the edge cells get longer). The target halffaces are coloured by
their area deviation relative to the start: blue is as far from the target as at the start, white
is at the target.

![Arearise chosen halffaces of a volmesh to a target area](../assets/images/examples/multi_cell_polyhedron_01_50_volmesh_arearisation.png)

```python title="examples/multi_cell_polyhedron/01_50_volmesh_arearisation.py"
--8<-- "examples/multi_cell_polyhedron/01_50_volmesh_arearisation.py"
```
