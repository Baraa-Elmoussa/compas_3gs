# API Reference

```python
from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer  # needs compas_viewer
from compas_3gs.diagrams import Cell, EGI, FormNetwork, ForceVolMesh, CellNetwork
from compas_3gs.algorithms import truss_equilibrium, cell_from_forces, volmesh_reciprocate
```

- [compas_3gs.model](compas_3gs.model.md): `GSModel`, a structure with supports and loads, its equilibrium and its force polyhedra.
- [compas_3gs.viewer](compas_3gs.viewer.md): `GSViewer`, the viewer of a model.
- [compas_3gs.diagrams](compas_3gs.diagrams.md): the form and force diagrams: cells, EGIs, form networks, polyhedral force diagrams.
- [compas_3gs.algorithms](compas_3gs.algorithms.md): EGI, cell construction, planarisation, arearisation, duality, reciprocation, unified diagrams, equilibrium.
- [compas_3gs.operations](compas_3gs.operations.md): topological and geometric operations on cells and force diagrams.
- [compas_3gs.datastructures](compas_3gs.datastructures.md): the mesh, network and volmesh base classes.
- [compas_3gs.utilities](compas_3gs.utilities.md): force magnitudes and colours, flatness and areaness, pairing of halffaces and edges.
