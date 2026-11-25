# Examples

Each example is a complete script in the `examples` folder of the repository, drawn with compas_viewer.
Run them from the root of the repository:

```bash
python examples/model/01_tripod.py
```

## Models

Structures with `GSModel` and `GSViewer`: supports, loads, equilibrium, force polyhedra.

<div class="grid cards" markdown>

-   [![A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view](../assets/images/examples/model_01_tripod.png)](model_01_tripod.md)

    **[A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view](model_01_tripod.md)**

-   [![A tetrahelix tower on a pinned, a slider and a roller support, with `GSModel`](../assets/images/examples/model_02_space_truss.png)](model_02_space_truss.md)

    **[A tetrahelix tower on a pinned, a slider and a roller support, with `GSModel`](model_02_space_truss.md)**

-   [![Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces](../assets/images/examples/model_03_force_first.png)](model_03_force_first.md)

    **[Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces](model_03_force_first.md)**

-   **[What `GSModel.solve` does when statics alone has many solutions, or none](model_04_determinacy.md)**

-   [![A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model](../assets/images/examples/model_05_warren_truss_fem.png)](model_05_warren_truss_fem.md)

    **[A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model](model_05_warren_truss_fem.md)**

-   [![A lattice transmission tower with two cross-arms, under conductor and wind loads](../assets/images/examples/model_06_transmission_tower.png)](model_06_transmission_tower.md)

    **[A lattice transmission tower with two cross-arms, under conductor and wind loads](model_06_transmission_tower.md)**

-   [![A Schwedler dome: rings, meridians and one diagonal per panel, on a ring of pinned supports](../assets/images/examples/model_07_schwedler_dome.png)](model_07_schwedler_dome.md)

    **[A Schwedler dome: rings, meridians and one diagonal per panel, on a ring of pinned supports](model_07_schwedler_dome.md)**

-   [![A tensegrity prism: three struts held apart by nine cables, in equilibrium without any load](../assets/images/examples/model_08_tensegrity_prism.png)](model_08_tensegrity_prism.md)

    **[A tensegrity prism: three struts held apart by nine cables, in equilibrium without any load](model_08_tensegrity_prism.md)**

</div>

## Physical structures

The same kind of structures built step by step with the algorithms: equilibrium, EGI, cells, duality and reciprocation.

<div class="grid cards" markdown>

-   [![A tripod: three legs carrying a load at the apex, and its force polyhedron](../assets/images/examples/physical_structures_10_tripod.png)](physical_structures_10_tripod.md)

    **[A tripod: three legs carrying a load at the apex, and its force polyhedron](physical_structures_10_tripod.md)**

-   [![A space truss tower (a tetrahelix) and one force polyhedron per node](../assets/images/examples/physical_structures_20_space_truss.png)](physical_structures_20_space_truss.md)

    **[A space truss tower (a tetrahelix) and one force polyhedron per node](physical_structures_20_space_truss.md)**

-   [![Force first: design a polyhedral force diagram, get the space frame it stands for](../assets/images/examples/physical_structures_30_force_first_frames.png)](physical_structures_30_force_first_frames.md)

    **[Force first: design a polyhedral force diagram, get the space frame it stands for](physical_structures_30_force_first_frames.md)**

-   [![Thrust line of a semicircular masonry arch, and what 3GS can and cannot do with it](../assets/images/examples/physical_structures_40_arch_thrust_line.png)](physical_structures_40_arch_thrust_line.md)

    **[Thrust line of a semicircular masonry arch, and what 3GS can and cannot do with it](physical_structures_40_arch_thrust_line.md)**

</div>

## Polyhedral cell

Constructing and manipulating a single polyhedral cell, the force polyhedron of one node.

<div class="grid cards" markdown>

-   [![Pull a face of a cell along its normal](../assets/images/examples/cell_manipulation_00_10_cell_pull_face.png)](cell_manipulation_00_10_cell_pull_face.md)

    **[Pull a face of a cell along its normal](cell_manipulation_00_10_cell_pull_face.md)**

-   [![Arearise one face of a cell: move it along its normal until it reaches a target area](../assets/images/examples/cell_manipulation_00_20_cell_arearise_face.png)](cell_manipulation_00_20_cell_arearise_face.md)

    **[Arearise one face of a cell: move it along its normal until it reaches a target area](cell_manipulation_00_20_cell_arearise_face.md)**

-   [![Construct a cell from force vectors in equilibrium](../assets/images/examples/cell_manipulation_00_30_cell_constructor.png)](cell_manipulation_00_30_cell_constructor.md)

    **[Construct a cell from force vectors in equilibrium](cell_manipulation_00_30_cell_constructor.md)**

-   [![Planarise the faces of a cell](../assets/images/examples/cell_manipulation_00_40_cell_planarisation.png)](cell_manipulation_00_40_cell_planarisation.md)

    **[Planarise the faces of a cell](cell_manipulation_00_40_cell_planarisation.md)**

-   [![Arearise several faces of a cell at once, optionally keeping all face normals fixed](../assets/images/examples/cell_manipulation_00_41_cell_arearisation.png)](cell_manipulation_00_41_cell_arearisation.md)

    **[Arearise several faces of a cell at once, optionally keeping all face normals fixed](cell_manipulation_00_41_cell_arearisation.md)**

</div>

## Multi-cell polyhedron

Polyhedral force diagrams of several cells, their dual form diagrams, reciprocation, planarisation and visualisation.

<div class="grid cards" markdown>

-   [![Construct a force volmesh from a set of cells and show its cells, vertices and vertex keys](../assets/images/examples/multi_cell_polyhedron_01_10_volmesh_constructor.png)](multi_cell_polyhedron_01_10_volmesh_constructor.md)

    **[Construct a force volmesh from a set of cells and show its cells, vertices and vertex keys](multi_cell_polyhedron_01_10_volmesh_constructor.md)**

-   [![Construct the dual volmesh of a force volmesh](../assets/images/examples/multi_cell_polyhedron_01_20_volmesh_duality_volmesh.png)](multi_cell_polyhedron_01_20_volmesh_duality_volmesh.md)

    **[Construct the dual volmesh of a force volmesh](multi_cell_polyhedron_01_20_volmesh_duality_volmesh.md)**

-   [![Make the dual network of a force volmesh and show how their elements correspond](../assets/images/examples/multi_cell_polyhedron_01_21_volmesh_duality_network.png)](multi_cell_polyhedron_01_21_volmesh_duality_network.md)

    **[Make the dual network of a force volmesh and show how their elements correspond](multi_cell_polyhedron_01_21_volmesh_duality_network.md)**

-   [![Reciprocate a polyhedral force diagram and its form diagram](../assets/images/examples/multi_cell_polyhedron_01_30_volmesh_reciprocation.png)](multi_cell_polyhedron_01_30_volmesh_reciprocation.md)

    **[Reciprocate a polyhedral force diagram and its form diagram](multi_cell_polyhedron_01_30_volmesh_reciprocation.md)**

-   [![Planarise the faces of a volmesh, with some of its vertices fixed](../assets/images/examples/multi_cell_polyhedron_01_40_volmesh_planarisation.png)](multi_cell_polyhedron_01_40_volmesh_planarisation.md)

    **[Planarise the faces of a volmesh, with some of its vertices fixed](multi_cell_polyhedron_01_40_volmesh_planarisation.md)**

-   [![Arearise chosen halffaces of a volmesh to a target area](../assets/images/examples/multi_cell_polyhedron_01_50_volmesh_arearisation.png)](multi_cell_polyhedron_01_50_volmesh_arearisation.md)

    **[Arearise chosen halffaces of a volmesh to a target area](multi_cell_polyhedron_01_50_volmesh_arearisation.md)**

-   [![Modify a polyhedral force diagram: set vertex fixities, move a vertex and lift a boundary vertex](../assets/images/examples/multi_cell_polyhedron_01_60_volmesh_modification.png)](multi_cell_polyhedron_01_60_volmesh_modification.md)

    **[Modify a polyhedral force diagram: set vertex fixities, move a vertex and lift a boundary vertex](multi_cell_polyhedron_01_60_volmesh_modification.md)**

-   [![Colour the corresponding elements of a polyhedral form and force diagram](../assets/images/examples/multi_cell_polyhedron_01_70_volmesh_color_id.png)](multi_cell_polyhedron_01_70_volmesh_color_id.md)

    **[Colour the corresponding elements of a polyhedral form and force diagram](multi_cell_polyhedron_01_70_volmesh_color_id.md)**

-   [![Five ways to visualise a polyhedral force diagram and its reciprocal form diagram](../assets/images/examples/multi_cell_polyhedron_01_71_volmesh_visualisation.png)](multi_cell_polyhedron_01_71_volmesh_visualisation.md)

    **[Five ways to visualise a polyhedral force diagram and its reciprocal form diagram](multi_cell_polyhedron_01_71_volmesh_visualisation.md)**

-   [![Unified diagram of a multi-cell force diagram and its form diagram](../assets/images/examples/multi_cell_polyhedron_01_80_volmesh_unified_diagram.png)](multi_cell_polyhedron_01_80_volmesh_unified_diagram.md)

    **[Unified diagram of a multi-cell force diagram and its form diagram](multi_cell_polyhedron_01_80_volmesh_unified_diagram.md)**

</div>
