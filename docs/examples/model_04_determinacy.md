# What `GSModel.solve` does when statics alone has many solutions, or none

- Four legs to pinned feet: one leg too many, the truss is statically indeterminate with one
  self-stress state. Every state in equilibrium is one of them plus any multiple of it; `solve` gives
  the elastic one for the member stiffness, or the one fixed by a prescribed force.

- Two legs: the apex can swing about the line through the feet, a mechanism under a sideways load.

- The same two legs carry a load in their own plane: a mechanism that carries this load only.

```python title="examples/model/04_determinacy.py"
--8<-- "examples/model/04_determinacy.py"
```
