# User Guide

The online documentation of `compas_3gs` is organised as follows.

<br>


------------------------------------------------------------------------

## Theoretical background

This section of the documentation provides a brief theoretical background of graphic statics.
The main motivations for the development of `compas_3gs` is also stated.
In addition, the intended design workflow using the package is described.
The basic nomenclature and terminology are also defined in this section.

- [Theoretical background](theory/index.md)

<br>


------------------------------------------------------------------------

## Models

These pages describe `GSModel` and `GSViewer`: a structure from its member lines, its supports and
loads, its equilibrium (determinate or not), the force polyhedra of its nodes, and the viewer.

- [The model](model.md): form first and force first, supports, loads, saving.
- [Solving](solving.md): statics, self-stress states, and how to pick a state of an indeterminate truss.
- [Force polyhedra](force_polyhedra.md): the force polyhedron of every node, and which nodes have one.
- [Viewing](viewing.md): the groups, the info panel, and changing the forces while viewing.

<br>


------------------------------------------------------------------------

## User manual

In this section, various examples are presented to demonstrate how the package can be used in design scenarios.
For each example, the underlying concepts of the computational methods used are visually expalined using diagrams.
The full source code for each example is provided, as well as the relevant CAD files.

- [Polyhedral cell](polyhedral_cell/index.md)
- [Multi-cell polyhedron](multi_cell_polyhedron/index.md)
- [Cell network](cell_network.md)
- [Materialisation](materialisation.md)

<br>


------------------------------------------------------------------------

## API Reference

This section explains in detail all the class objects, algorithms and functionalities of the package.

- [API Reference](../api/index.md)
