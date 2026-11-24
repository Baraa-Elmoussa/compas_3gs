# TODO

- Upstream compas: `VolMesh.cell_vertex_faces` / `cell_vertex_neighbors` index `dict_keys` and raise `TypeError` in 2.15.1; fix with `list(...)`. compas_3gs goes through `cell_to_mesh` meanwhile (`EGI.from_volmesh_cell`). Groups 3, 10.
- `volmesh_cell_split_vertex` reads `self.c_data` (compas 0.x) and has never run under compas 1 or 2; decide: port with `EGI.from_volmesh_cell` + `split_face`, or raise `NotImplementedError`. Group 10.
- `cell_cull_zero_faces` (`operations/celloperations.py`) is an empty stub, not in `__all__`, with no callers; decide: implement or delete.
