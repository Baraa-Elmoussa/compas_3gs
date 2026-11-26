# AGENTS.md

`compas_3gs` is 3D graphic statics for COMPAS 2: the equilibrium of a structure drawn as two
reciprocal diagrams. The **form diagram** is the structure (a network of members). The **force
diagram** is a set of polyhedral cells, one per node, with one face per force on that node: the face
is perpendicular to the force and its area is the force's magnitude.

This file is the entry point for agents. The package already implements the statics and the
geometry of 3D graphic statics; use its API instead of writing your own equilibrium solver, EGI,
planarisation or reciprocation. The files in [`.agents/`](.agents/) list that API by task.

## The theory

Read these pages of the documentation (Markdown, in `docs/`) before changing an algorithm:

- [Graphic statics](docs/userguide/theory/00_graphic_statics.md), [Motivation](docs/userguide/theory/01_motivation.md)
- [Design workflow](docs/userguide/theory/02_design_workflow.md), [Implementation](docs/userguide/theory/03_implementation.md)
- [Nomenclature](docs/userguide/theory/04_nomenclature.md): the names used in the code (cell, EGI, halfface, ...)
- The User Guide, one page per method: [polyhedral cell](docs/userguide/polyhedral_cell/index.md),
  [multi-cell polyhedron](docs/userguide/multi_cell_polyhedron/index.md),
  [cell network](docs/userguide/cell_network.md), [materialisation](docs/userguide/materialisation.md)

## What the package can do

- **Model a truss** from member lines, with supports and forces, and solve its equilibrium:
  determinate trusses by statics; indeterminate ones by member stiffness (the elastic state) or by
  prescribed forces or reactions (any state in equilibrium); mechanisms are detected.
  → [`.agents/model.md`](.agents/model.md)
- **Build force polyhedra**: one cell per node from its forces (form first), or a polyhedral force
  diagram whose dual network is the form diagram (force first).
  → [`.agents/model.md`](.agents/model.md), [`.agents/algorithms.md`](.agents/algorithms.md)
- **Construct and transform cells and multi-cell force diagrams**: EGI, planarisation to target
  normals and areas, arearisation, duality, reciprocation, unified diagrams, topological operations.
  → [`.agents/algorithms.md`](.agents/algorithms.md), [`.agents/operations.md`](.agents/operations.md),
  [`.agents/diagrams.md`](.agents/diagrams.md)
- **Draw** a model in compas_viewer, grouped, with compression and tension, supports, loads,
  reactions, force cells and an info panel with the force of a selected member and the forces on
  the node of a selected load, support, reaction or cell.
  → [`.agents/viewer.md`](.agents/viewer.md)

[`.agents/tasks.md`](.agents/tasks.md) maps common tasks to the call and to an example that does it.

## How to work here

1. Look up the task in [`.agents/tasks.md`](.agents/tasks.md), then the function in the topic file.
2. Start from the closest script in `examples/` (`examples/model/` for the `GSModel` API); every
   example runs from the repository root and is tested by `tests/test_examples.py`.
3. Read [`.agents/compas2.md`](.agents/compas2.md) before touching a `VolMesh`: several compas
   2.15 and 2.16 methods are broken or surprising, and the package works around them.
4. Conventions: member forces are positive in tension and negative in compression; no units are
   imposed (the examples use kN and m); node and edge keys are not indices.

Setup and checks (Python 3.12, `compas >=2.15,<3`; `compas_viewer` only for drawing):

```bash
pip install -e ".[dev]"
pytest
ruff check . && ruff format --check .
mkdocs build --strict
```

Open design questions and known upstream bugs are in [`TODO.md`](TODO.md); changes go in
[`CHANGELOG.md`](CHANGELOG.md).
