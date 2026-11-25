# Construct a cell from force vectors in equilibrium

The vectors give an EGI, the EGI a unit cell, and the cell is planarised to face areas equal
to the vector lengths and face normals equal to the vector directions. Faces that come from
crossings of the EGI arcs (zero faces) get a target area of 0.

![Construct a cell from force vectors in equilibrium](../assets/images/examples/cell_manipulation_00_30_cell_constructor.png)

```python title="examples/cell_manipulation/00_30_cell_constructor.py"
--8<-- "examples/cell_manipulation/00_30_cell_constructor.py"
```
