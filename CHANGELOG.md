# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

* Added the sample data in `data/` as plain JSON (cells as vertices and faces, volmeshes as vertices and cells), converted from the `.3dm` files by `scripts/convert_3dm.py`.
* Added examples for compas_viewer: `examples/cell_manipulation/` (5), `examples/multi_cell_polyhedron/` (10) and `examples/physical_structures/` (4).
* Added type hints and NumPy-style Markdown docstrings to the ported classes and functions.
* Added `TODO.md` with the open design questions and the known compas 2.15.1 and compas_viewer 2.0.2 issues.
* Added `pyproject.toml` and `tasks.py` with `compas_invocations2`.
* Added supports, applied forces and member forces to `compas_3gs.diagrams.FormNetwork`: `node_at`, `set_support`, `add_force`, `supports`, `node_load`, `node_reaction`, `edge_force` and `node_forces`, stored in the attributes `support`, `px`, `py`, `pz`, `_rx`, `_ry`, `_rz` and `_f`.
* Added `compas_3gs.algorithms.truss_equilibrium` and `EquilibriumResult`: member forces and support reactions of a pin-jointed truss from equilibrium.
* Added `compas_3gs.algorithms.cell_from_forces`: the force polyhedron of a node from the forces on it.
* Added `compas_3gs.model.GSModel`: a form diagram with supports and forces, its equilibrium (`solve`) and the force polyhedron of every node (`derive_force_diagram`), built from member lines (`from_lines`) or from a polyhedral force diagram (`from_force_diagram`).
* Added `compas_3gs.viewer.GSViewer`, a compas_viewer `Viewer` that draws a `GSModel` in groups: members in compression and tension with their forces, support glyphs, loads, reactions and the force diagram. The module needs compas_viewer.
* Added examples of `GSModel` and `GSViewer` in `examples/model/`.
* Added tests for every module in `tests/compas_3gs/`, mirroring `src/`, with JSON round trips of the data structures, and `tests/test_examples.py`, which runs every example (marked `gui`).

### Changed

* Updated workflow.
* Changed the documentation from Sphinx to MkDocs (Material and mkdocstrings), with the pages converted to Markdown unchanged and an API page per subpackage.
* Changed the requirements to compas 2.15 (`compas >=2.15, <3`) and Python 3.12.
* Changed `Mesh3gs`, `Cell` and `EGI` to compas 2 `Mesh`, `Network3gs`, `FormNetwork` and `CellNetwork` to compas 2 `Graph`, and `VolMesh3gs` and `ForceVolMesh` to compas 2 `VolMesh`, all with constructors of the form `__init__(self, *args, name=..., **kwargs)`.
* Changed `Mesh3gs.face_area` to the vector area, and `Mesh3gs.face_normal` to return a `Vector`, with a fallback for faces whose normal is zero.
* Changed `Cell.cell_center` to the volume centroid.
* Changed the EGI construction to work without Rhino (great-circle intersections, points compared by distance with `TOL.angular`).
* Changed the colour utilities to return 0–255 tuples from `compas.colors.Color`, with values clamped to 0–1.
* Changed `point_reflection` to work for graphs, meshes and volmeshes.
* Changed the volmesh operations to rebuild cells through the public API, keeping cell, face and edge attributes.
* Changed mutable default arguments to `None`.
* Changed the volmesh `scale` to `attributes["scale"]`.
* Moved the old Rhino examples and `.3dm` files to `examples/rhino/`.
* Fixed `add_vertex(vkey=...)` calls in `volmesh_dual_volmesh` and `volmesh_vertex_merge`, which stored the key as an attribute and let compas pick another key.
* Fixed the `CellNetwork` default attributes, which made all nodes share one `Cell`.
* Fixed `cell_arearise_face`, which called the undefined `face_oriented_area` and `face_oriented_normal`.
* Fixed `cell_center` for cells with gaps in their vertex keys.
* Fixed `volmesh_merge_adjacent_halffaces`, which could build an inward cell and rejected halffaces of neighbouring cells.
* Fixed the faces of `volmesh_dual_volmesh`, which were not in cyclic order around their edge.
* Fixed duplicate edges in `volmesh_dual_network`.
* Fixed `volmesh_planarise` targets keyed by the opposite halfface, which raised `KeyError` (normals, centres) or were ignored (areas).
* Fixed `volmesh_vertex_merge`, which kept a repeated vertex where the first and last vertex of a face were merged.
* Fixed `Network3gs.edge_avg_length`.
* Fixed `resultant_vector`, which ignored `locations` and placed the resultant at the weighted centroid of the vectors themselves.

### Removed

* Removed the Rhino plugin installer (`install.py`, `activate.py`, `scripts/pack.py`, `scripts/unzip_electron.py`) and the dependency on `compas_skeleton`.
* Removed `setup.py`, `setup.cfg`, `MANIFEST.in`, `.bumpversion.cfg` and `pytest.ini`.
* Removed the Rhino drawing methods of `Mesh3gs`.
* Removed the `Network3gs` aliases `vertices`, `vertex_coordinates` and `vertex_neighbors`, `Network3gs.bounding_box` and the `Network3gs.edge_vector` override.
* Removed `VolMesh3gs.delete_cell` (compas has it), `VolMesh3gs.cell_vertex_delete` and `VolMesh3gs.delete_halfface`.
* Removed `VolMesh3gs.boundary_halfface_manifold_neighbors` and `boundary_halfface_manifold_neighborhood`, which skipped neighbours in the same cell and passed the wrong key in their ring loop.
* Removed private helpers that nothing uses any more from `algorithms/planarisation.py` and `utilities/geometry.py`.
