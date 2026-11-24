"""Convert the Rhino example files in `examples/rhino/rhino_files` to COMPAS JSON in `data`.

The examples read only the JSON, so they run without Rhino. Run this script once, with
`rhino3dm` installed, whenever a `.3dm` file changes:

    python scripts/convert_3dm.py

Output (coordinates in the units of the `.3dm` file):

- `cell_gss.json`: `{"vertices": [[x, y, z], ...], "faces": [[i, j, k, ...], ...]}`, one closed
  polyhedral cell. Build it with `Cell.from_vertices_and_faces(data["vertices"], data["faces"])`.
- `cells.json`: `{"cells": [{"vertices": ..., "faces": ...}, ...], "point": [x, y, z]}`.
- `equilibrated_force_vectors.json`: a list of `{"origin": [x, y, z], "vectors": [[x, y, z], ...]}`,
  one per group of force vectors in equilibrium (lines that share a start point).
- `volmesh_*.json`: `{"vertices": [[x, y, z], ...], "cells": [[face, ...], ...]}`, one cell per
  polysurface, each face a list of vertex indices, with the vertices shared between cells. Build
  the diagram with `ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])`.

All files hold plain vertices and faces rather than serialised COMPAS objects: the files stay the
same when the script runs again (a serialised object gets a new `guid` every time), they do not
depend on the class names of compas_3gs, and in compas 2.15.1 a `VolMesh` read back from JSON (or
copied) gets every face once per edge.

Polysurface faces are read from the outer loop of each Brep face, not from the render mesh, so a
quad stays a quad. Corners of neighbouring cells are welded when they are closer than 1% of the
shortest edge in the file, or the model tolerance of the file if that is larger: the cells that
3GS baked to Rhino (`volmesh_donut`, `volmesh_honeycomb`) have shared corners up to 0.02 apart,
far more than the model tolerance, while their shortest edges are about 2 and 5 long.

"""

import os

import rhino3dm

import compas
from compas.datastructures import Mesh
from compas.datastructures import VolMesh
from compas.geometry import distance_point_point

HERE = os.path.dirname(__file__)
RHINO_FILES = os.path.abspath(os.path.join(HERE, "..", "examples", "rhino", "rhino_files"))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))

POLYSURFACE_FILES = [
    "cell_gss",
    "cells",
    "volmesh_cubes_4x4",
    "volmesh_donut",
    "volmesh_honeycomb",
    "volmesh_mat",
    "volmesh_shells",
    "volmesh_warped",
]

os.makedirs(DATA, exist_ok=True)

# ==============================================================================
# Polysurfaces: cells and volmeshes
# ==============================================================================

for stem in POLYSURFACE_FILES:
    model = rhino3dm.File3dm.Read(os.path.join(RHINO_FILES, stem + ".3dm"))

    objects = list(model.Objects)
    if all((obj.Attributes.Name or "").startswith("cell.") for obj in objects):
        objects.sort(key=lambda obj: int(obj.Attributes.Name.split(".")[1]))

    breps = []
    points = []
    for obj in objects:
        geometry = obj.Geometry
        if isinstance(geometry, rhino3dm.Point):
            points.append([geometry.Location.X, geometry.Location.Y, geometry.Location.Z])
        if isinstance(geometry, rhino3dm.Extrusion):
            geometry = geometry.ToBrep(True)
        if isinstance(geometry, rhino3dm.Brep):
            breps.append(geometry)

    shortest = min(
        brep.Vertices[trim.StartVertexIndex].Location.DistanceTo(brep.Vertices[trim.EndVertexIndex].Location)
        for brep in breps
        for face in brep.Faces
        for trim in face.OuterLoop.Trims
        if trim.StartVertexIndex != trim.EndVertexIndex
    )
    tol = max(model.Settings.ModelAbsoluteTolerance, 0.01 * shortest)

    vertices = []
    cells = []

    for geometry in breps:
        # weld the Brep vertices into the shared vertex list
        index = {}
        for vertex in geometry.Vertices:
            xyz = [vertex.Location.X, vertex.Location.Y, vertex.Location.Z]
            for i, other in enumerate(vertices):
                if distance_point_point(xyz, other) < tol:
                    index[vertex.VertexIndex] = i
                    break
            else:
                index[vertex.VertexIndex] = len(vertices)
                vertices.append(xyz)

        # one polygon per Brep face, from its outer loop
        faces = []
        for face in geometry.Faces:
            cycle = [index[trim.StartVertexIndex] for trim in face.OuterLoop.Trims]
            cycle = [v for i, v in enumerate(cycle) if v != cycle[i - 1]]
            if face.OrientationIsReversed:
                cycle.reverse()
            faces.append(cycle)

        # Rhino solids are outward, but check: a cell with negative volume is turned around
        used = sorted(set(v for face in faces for v in face))
        local = {v: i for i, v in enumerate(used)}
        mesh = Mesh.from_vertices_and_faces([vertices[v] for v in used], [[local[v] for v in face] for face in faces])
        if mesh.volume(unify_cycles=False) < 0:
            faces = [face[::-1] for face in faces]

        cells.append(faces)

    if stem.startswith("volmesh_"):
        volmesh = VolMesh.from_vertices_and_cells(vertices, cells)
        interior = len(list(volmesh.halffaces())) - len(volmesh.halffaces_on_boundaries())
        compas.json_dump({"vertices": vertices, "cells": cells}, os.path.join(DATA, stem + ".json"), pretty=True)
        print(f"{stem}: {len(cells)} cells, {len(vertices)} vertices, {interior} interior halffaces, weld tolerance {tol:.3g}")
        continue

    polyhedra = []
    for faces in cells:
        used = sorted(set(v for face in faces for v in face))
        local = {v: i for i, v in enumerate(used)}
        polyhedra.append({"vertices": [vertices[v] for v in used], "faces": [[local[v] for v in face] for face in faces]})

    if stem == "cell_gss":
        compas.json_dump(polyhedra[0], os.path.join(DATA, stem + ".json"), pretty=True)
        print(f"{stem}: 1 cell, {len(polyhedra[0]['vertices'])} vertices, {len(polyhedra[0]['faces'])} faces")
    else:
        compas.json_dump({"cells": polyhedra, "point": points[0] if points else None}, os.path.join(DATA, stem + ".json"), pretty=True)
        print(f"{stem}: {len(polyhedra)} cells, point {points[0] if points else None}")

# ==============================================================================
# Force vectors in equilibrium
# ==============================================================================

model = rhino3dm.File3dm.Read(os.path.join(RHINO_FILES, "equilibrated_force_vectors.3dm"))
tol = model.Settings.ModelAbsoluteTolerance

groups = []
for obj in model.Objects:
    line = obj.Geometry
    if not isinstance(line, rhino3dm.LineCurve):
        continue
    start = [line.PointAtStart.X, line.PointAtStart.Y, line.PointAtStart.Z]
    end = [line.PointAtEnd.X, line.PointAtEnd.Y, line.PointAtEnd.Z]
    vector = [end[0] - start[0], end[1] - start[1], end[2] - start[2]]
    for group in groups:
        if distance_point_point(group["origin"], start) < tol:
            group["vectors"].append(vector)
            break
    else:
        groups.append({"origin": start, "vectors": [vector]})

groups.sort(key=lambda group: group["origin"][1])
compas.json_dump(groups, os.path.join(DATA, "equilibrated_force_vectors.json"), pretty=True)
print(f"equilibrated_force_vectors: {len(groups)} groups, {[len(group['vectors']) for group in groups]} vectors")
