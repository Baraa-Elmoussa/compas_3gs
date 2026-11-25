# Arearise several faces of a cell at once, optionally keeping all face normals fixed

With fixed normals the targets must be compatible with a closed cell (the area-weighted normals
sum to zero); the targets below are reached with free normals.

The faces are coloured by their flatness relative to the start: blue is as warped as at the
start, white is flat.

![Arearise several faces of a cell at once, optionally keeping all face normals fixed](../assets/images/examples/cell_manipulation_00_41_cell_arearisation.png)

```python title="examples/cell_manipulation/00_41_cell_arearisation.py"
--8<-- "examples/cell_manipulation/00_41_cell_arearisation.py"
```
