# A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view

Form diagram: three legs from an apex to pinned feet, and a load at the apex.
`solve` gives the leg forces and the reactions from equilibrium; `derive_force_diagram` builds the
force polyhedron of every node that has at least four forces in three dimensions (here the apex:
three legs and the load). The feet have only a leg and a reaction, so they are skipped.

![A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view](../assets/images/examples/model_01_tripod.png)

```python title="examples/model/01_tripod.py"
--8<-- "examples/model/01_tripod.py"
```
