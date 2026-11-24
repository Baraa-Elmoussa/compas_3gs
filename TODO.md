# TODO

- Upstream compas: `VolMesh.cell_vertex_faces` / `cell_vertex_neighbors` index `dict_keys` and raise `TypeError` in 2.15.1; fix with `list(...)`. compas_3gs goes through `cell_to_mesh` meanwhile (`EGI.from_volmesh_cell`). Groups 3, 9.
- `volmesh_cell_split_vertex` reads `self.c_data` (compas 0.x) and has never run under compas 1 or 2; decide: port with `EGI.from_volmesh_cell` + `split_face`, or raise `NotImplementedError`. Group 9.
- Upstream compas: `VolMesh.__data__` writes each face of a cell once per halfedge (2.15.1), so a `VolMesh` read back from JSON or `copy()` gets every face several times (2×2×2 meshgrid: 192 halffaces instead of 48); It also raises `TypeError` as soon as any cell attribute is set, because its `"cell_data"` line unpacks the dict keys as pairs. compas_3gs never copies or serialises a volmesh meanwhile: `data/volmesh_*.json` is plain vertices and cells. Groups 7, 9, 14.
- `cell_cull_zero_faces` (`operations/celloperations.py`) is an empty stub, not in `__all__`, with no callers; decide: implement or delete.
- `docs/user_manual` links point at the old `examples/rhino_files` and `examples/_old_examples` paths (now `examples/rhino/rhino_files` and `examples/rhino/old_examples`); fix them in Group 19.
