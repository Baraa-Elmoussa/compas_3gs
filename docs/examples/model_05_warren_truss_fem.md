# A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model

The truss has a triangular section: two bottom chords, one top chord, Warren diagonals from each
top node to the four bottom nodes of its bay, cross ties and lateral bracing in the bottom plane.
It carries a deck load on the bottom nodes, a wind load on the top nodes and a vehicle on one side
of the deck.

A. On three bearings (pinned, roller, slider) the truss is statically determinate: as many
   independent equilibrium equations as unknowns. 3GS gets the member forces from statics alone,
   so the finite element model must give the same forces for any member stiffness. It also gives
   what statics cannot: the displacements. 3GS also gives what the finite element model does not:
   the force polyhedron of every node, here kept only where it reaches the forces (`check`).

B. With a fourth bearing the truss is statically indeterminate, with one self-stress state: every
   state in equilibrium with the loads is one of them plus any multiple of it, a line of states.

   - The member stiffness picks one point on that line, the elastic state. 3GS gives it from
     statics and stiffness (the least complementary energy), and it is the finite element state.
     Here the bracing of the first half of the span is `SOFT_BRACING` times softer.

   - Any other point is a state in equilibrium too: pick it in the input with `FOURTH_BEARING`, the
     vertical reaction at the fourth bearing, or move the slider in the viewer. The script prints
     the range in which no bearing lifts off, and builds the force polyhedra of the chosen state.
     At 0 kN the fourth bearing carries nothing, and the state is that of A.
   (Scaling a whole member group, or any change under loads symmetric about the middle of the span,
   would not move the elastic state here: the self-stress state twists the truss, with the bearing
   reactions alternating around the corners.)

The finite element model is the direct stiffness method for pin-jointed bars: each member adds
EA/L * d d^T to the stiffness matrix (d its unit direction), the free degrees of freedom are solved
for the loads, and the axial force of a member is EA/L times its elongation.

![A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model](../assets/images/examples/model_05_warren_truss_fem.png)

```python title="examples/model/05_warren_truss_fem.py"
--8<-- "examples/model/05_warren_truss_fem.py"
```
