# A tripod: three legs carrying a load at the apex, and its force polyhedron

Form diagram: the apex node, three legs to pinned feet, and an applied load.
Force diagram: one closed polyhedral cell for the apex node. Each face belongs to one force at
the node: its normal is the direction of the force, its area is the magnitude. Because the forces
are in equilibrium, the cell closes (the area-weighted face normals sum to zero).

Steps:
1. Solve the apex equilibrium for the three leg forces (3 equations, 3 unknowns).
2. Give the four force vectors (3 legs + load) to 3GS: EGI -> unit cell -> planarise to areas.
3. Check: every face area equals its force, every face normal is parallel to its force (so each
   face is perpendicular to its member: the reciprocal relation of graphic statics).

![A tripod: three legs carrying a load at the apex, and its force polyhedron](../assets/images/examples/physical_structures_10_tripod.png)

```python title="examples/physical_structures/10_tripod.py"
--8<-- "examples/physical_structures/10_tripod.py"
```
