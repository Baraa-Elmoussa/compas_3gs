import json
import os
import random

import pytest

from compas.geometry import Box
from compas.geometry import add_vectors
from compas.geometry import distance_point_plane
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs import DATA
from compas_3gs.algorithms import cell_from_egi
from compas_3gs.algorithms import cell_planarise
from compas_3gs.algorithms import egi_from_vectors
from compas_3gs.algorithms import volmesh_planarise
from compas_3gs.diagrams import Cell
from compas_3gs.diagrams import ForceVolMesh

# the default convergence tolerances of cell_planarise and volmesh_planarise
TOLERANCE_FLAT = 0.001
TOLERANCE_AREA = 0.001
TOLERANCE_PERP = 0.001

with open(os.path.join(DATA, "cells.json")) as f:
    CELLS = json.load(f)["cells"]

with open(os.path.join(DATA, "equilibrated_force_vectors.json")) as f:
    FORCE_GROUPS = json.load(f)


@pytest.fixture
def warped_cube(request):
    """A unit cube whose vertices are moved by up to 0.1 along each axis, with the seed `request.param`."""
    cell = Cell.from_shape(Box(1))
    rng = random.Random(request.param)
    for vkey in cell.vertices():
        cell.vertex_attributes(vkey, "xyz", add_vectors(cell.vertex_coordinates(vkey), [rng.uniform(-0.1, 0.1) for _ in range(3)]))
    return cell


@pytest.fixture
def warped_meshgrid(request):
    """A 2x2x2 meshgrid of unit cubes whose vertices are moved by up to 0.2 along each axis, with the seed `request.param`."""
    volmesh = ForceVolMesh.from_meshgrid(dx=2, nx=2)
    rng = random.Random(request.param)
    for vkey in volmesh.vertices():
        volmesh.vertex_attributes(vkey, "xyz", add_vectors(volmesh.vertex_coordinates(vkey), [rng.uniform(-0.2, 0.2) for _ in range(3)]))
    return volmesh


# ==============================================================================
# cell_planarise
# ==============================================================================


@pytest.mark.parametrize("warped_cube", [0, 1, 2], indirect=True)
def test_cell_planarise_warped_cube(warped_cube):
    cell = warped_cube
    flatness = [distance_point_plane(cell.vertex_coordinates(v), (cell.face_centroid(f), cell.face_normal(f))) for f in cell.faces() for v in cell.face_vertices(f)]
    assert max(flatness) > TOLERANCE_FLAT

    calls = []
    cell_planarise(cell, kmax=100, callback=lambda cell, k, args: calls.append((k, args)), callback_args="args")

    # the callback runs after every iteration but the one that converges
    assert len(calls) < 100
    assert calls == [(k, "args") for k in range(len(calls))]
    assert cell.number_of_vertices() == 8
    assert cell.number_of_edges() == 12
    assert cell.number_of_faces() == 6
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    center = cell.centroid()
    for fkey in cell.faces():
        normal = cell.face_normal(fkey)
        centroid = cell.face_centroid(fkey)
        assert dot_vectors(normal, subtract_vectors(centroid, center)) > 0
        for vkey in cell.face_vertices(fkey):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (centroid, normal)) < TOLERANCE_FLAT


def test_cell_planarise_target_areas():
    cell = Cell.from_vertices_and_faces(CELLS[0]["vertices"], CELLS[0]["faces"])
    target_areas = {2: 50.0, 3: 70.0, 6: 15.0}

    calls = []
    cell_planarise(cell, kmax=500, target_areas=target_areas, callback=lambda cell, k, args: calls.append(k))

    assert len(calls) < 500
    for fkey, target in target_areas.items():
        assert abs(cell.face_area(fkey) - target) < TOLERANCE_AREA
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    for fkey in cell.faces():
        for vkey in cell.face_vertices(fkey):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (cell.face_centroid(fkey), cell.face_normal(fkey))) < TOLERANCE_FLAT


def test_cell_planarise_cube_from_egi():
    axes = {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1], 5: [0, 0, -1]}
    # a cube of side 2/3
    cell = cell_from_egi(egi_from_vectors(axes, [0, 0, 0]))

    calls = []
    cell_planarise(cell, kmax=100, target_areas={fkey: 4.0 for fkey in axes}, target_normals=axes, callback=lambda cell, k, args: calls.append(k))

    # a cube of side 2
    assert len(calls) < 100
    for fkey in cell.faces():
        assert abs(cell.face_area(fkey) - 4.0) < TOLERANCE_AREA
        assert 1 - dot_vectors(cell.face_normal(fkey), axes[fkey]) < TOL.approximation
    assert volume_polyhedron(cell.to_vertices_and_faces()) == pytest.approx(8.0, rel=TOL.approximation)


@pytest.mark.parametrize("index", range(len(FORCE_GROUPS)))
def test_cell_planarise_equilibrated_force_vectors(index):
    # the faces of the cell get the areas and directions of the vectors; the faces from arc crossings (zero faces) get area 0
    origin = FORCE_GROUPS[index]["origin"]
    vectors = dict(enumerate(FORCE_GROUPS[index]["vectors"]))
    egi = egi_from_vectors(vectors, origin)
    cell = cell_from_egi(egi)
    target_areas = {fkey: length_vector(vectors[fkey]) if fkey in vectors else 0.0 for fkey in cell.faces()}
    target_normals = {fkey: egi.vertex_attribute(fkey, "normal") for fkey in cell.faces()}

    calls = []
    cell_planarise(cell, kmax=1000, target_areas=target_areas, target_normals=target_normals, collapse_edge_length=0.1, callback=lambda cell, k, args: calls.append(k))

    assert len(calls) < 1000
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    for fkey in cell.faces():
        assert abs(cell.face_area(fkey) - target_areas[fkey]) < TOLERANCE_AREA
    for fkey in vectors:
        assert 1 - dot_vectors(cell.face_normal(fkey), normalize_vector(vectors[fkey])) < TOL.approximation
        for vkey in cell.face_vertices(fkey):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (cell.face_centroid(fkey), cell.face_normal(fkey))) < TOLERANCE_FLAT


@pytest.mark.parametrize("warped_cube", [3], indirect=True)
def test_cell_planarise_fixed_vertices(warped_cube):
    cell = warped_cube
    cell.vertex_attribute(3, "z_fix", True)
    fixed = {vkey: cell.vertex_coordinates(vkey) for vkey in (0, 6)}
    xyz = cell.vertex_coordinates(3)

    cell_planarise(cell, kmax=100, fix_vkeys=list(fixed))

    for vkey in fixed:
        assert cell.vertex_coordinates(vkey) == fixed[vkey]
    assert cell.vertex_attribute(3, "z") == xyz[2]
    assert cell.vertex_coordinates(3)[:2] != xyz[:2]
    for fkey in cell.faces():
        for vkey in cell.face_vertices(fkey):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (cell.face_centroid(fkey), cell.face_normal(fkey))) < TOLERANCE_FLAT


def test_cell_planarise_avg_fkeys():
    cell = Cell.from_shape(Box(1))
    cell.vertex_attributes(0, "xyz", add_vectors(cell.vertex_coordinates(0), [0.1, 0.05, -0.1]))
    avg_fkeys = cell.vertex_faces(0)
    assert max(cell.face_area(fkey) for fkey in avg_fkeys) - min(cell.face_area(fkey) for fkey in avg_fkeys) > 0.05

    cell_planarise(cell, kmax=100, avg_fkeys=avg_fkeys)

    mean = sum(cell.face_area(fkey) for fkey in avg_fkeys) / len(avg_fkeys)
    for fkey in avg_fkeys:
        assert abs(cell.face_area(fkey) - mean) < TOLERANCE_AREA


def test_cell_planarise_callback_not_callable():
    with pytest.raises(Exception, match="not callable"):
        cell_planarise(Cell.from_shape(Box(1)), callback=1)


# ==============================================================================
# volmesh_planarise
# ==============================================================================


@pytest.mark.parametrize("warped_meshgrid", [0, 1], indirect=True)
def test_volmesh_planarise_warped_meshgrid(warped_meshgrid):
    volmesh = warped_meshgrid
    interior = [hfkey for hfkey in volmesh.faces() if not volmesh.is_halfface_on_boundary(hfkey)]
    assert len(interior) == 12
    for hfkey in interior:
        plane = (volmesh.halfface_center(hfkey), volmesh.halfface_normal(hfkey))
        assert max(distance_point_plane(volmesh.vertex_coordinates(vkey), plane) for vkey in volmesh.halfface_vertices(hfkey)) > TOLERANCE_FLAT

    calls = []
    volmesh_planarise(volmesh, kmax=200, callback=lambda volmesh, k, args: calls.append((k, args)), callback_args="args")

    assert len(calls) < 200
    assert calls == [(k, "args") for k in range(len(calls))]
    # every face, interior faces included, and both of its halffaces
    for hfkey in volmesh.halffaces():
        for vkey in volmesh.halfface_vertices(hfkey):
            assert distance_point_plane(volmesh.vertex_coordinates(vkey), (volmesh.halfface_center(hfkey), volmesh.halfface_normal(hfkey))) < TOLERANCE_FLAT
    # valid: closed cells, halffaces pointing out of their cell (signed volume), interior halffaces paired
    assert volmesh.number_of_cells() == 8
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0
    # the halffaces with the same vertices are each other's opposite; the others are the 24 boundary halffaces
    groups = {}
    for hfkey in volmesh.halffaces():
        groups.setdefault(frozenset(volmesh.halfface_vertices(hfkey)), []).append(hfkey)
    assert sorted(len(group) for group in groups.values()) == [1] * 24 + [2] * 12
    assert sorted(group[0] for group in groups.values() if len(group) == 1) == sorted(volmesh.halffaces_on_boundaries())
    for group in groups.values():
        if len(group) == 2:
            assert volmesh.halfface_opposite_halfface(group[0]) == group[1]
            assert volmesh.halfface_opposite_halfface(group[1]) == group[0]


def test_volmesh_planarise_warped_interior_face():
    volmesh = ForceVolMesh.from_meshgrid(dx=2, nx=2)
    interior_vkeys = [vkey for vkey in volmesh.vertices() if not volmesh.is_vertex_on_boundary(vkey)]
    assert len(interior_vkeys) == 1
    volmesh.vertex_attributes(interior_vkeys[0], "xyz", add_vectors(volmesh.vertex_coordinates(interior_vkeys[0]), [0.3, 0.2, -0.25]))
    warped = []
    for hfkey in volmesh.faces():
        plane = (volmesh.halfface_center(hfkey), volmesh.halfface_normal(hfkey))
        if max(distance_point_plane(volmesh.vertex_coordinates(vkey), plane) for vkey in volmesh.halfface_vertices(hfkey)) > TOLERANCE_FLAT:
            warped.append(hfkey)
    assert len(warped) == 12
    assert not any(volmesh.is_halfface_on_boundary(hfkey) for hfkey in warped)

    calls = []
    volmesh_planarise(volmesh, kmax=100, callback=lambda volmesh, k, args: calls.append(k))

    assert len(calls) < 100
    for hfkey in volmesh.faces():
        for vkey in volmesh.halfface_vertices(hfkey):
            assert distance_point_plane(volmesh.vertex_coordinates(vkey), (volmesh.halfface_center(hfkey), volmesh.halfface_normal(hfkey))) < TOLERANCE_FLAT
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0


@pytest.mark.parametrize("keyed_by", ["halfface", "opposite"])
def test_volmesh_planarise_targets_keyed_by_either_halfface(keyed_by):
    volmesh = ForceVolMesh.from_meshgrid(dx=2, nx=2)
    hfkey = next(hfkey for hfkey in volmesh.faces() if not volmesh.is_halfface_on_boundary(hfkey))
    opposite = volmesh.halfface_opposite_halfface(hfkey)
    # `faces()` yields one halfface per face: the opposite one is only reached through the targets
    assert opposite not in set(volmesh.faces())
    key = hfkey if keyed_by == "halfface" else opposite
    # one plane, tilted and moved from the face, and a smaller area
    normal = normalize_vector(add_vectors(volmesh.halfface_normal(hfkey), [0.2, 0.1, 0.15]))
    center = add_vectors(volmesh.halfface_center(hfkey), scale_vector(volmesh.halfface_normal(hfkey), 0.1))
    area = 0.8

    calls = []
    volmesh_planarise(volmesh, kmax=500, target_centers={key: center}, target_normals={key: normal}, target_areas={key: area}, callback=lambda volmesh, k, args: calls.append(k))

    assert len(calls) < 500
    assert 1 - abs(dot_vectors(volmesh.halfface_normal(hfkey), normal)) < TOLERANCE_PERP
    assert abs(volmesh.halfface_area(hfkey) - area) < TOLERANCE_AREA
    assert TOL.is_close(volmesh.halfface_area(opposite), volmesh.halfface_area(hfkey))
    for vkey in volmesh.halfface_vertices(hfkey):
        assert distance_point_plane(volmesh.vertex_coordinates(vkey), (center, normal)) < TOLERANCE_FLAT
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0


def test_volmesh_planarise_fixed_vertices():
    with open(os.path.join(DATA, "volmesh_warped.json")) as f:
        data = json.load(f)
    volmesh = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])
    # the six lowest vertices, as if the diagram stood on them
    zmin = min(volmesh.vertex_attribute(vkey, "z") for vkey in volmesh.vertices())
    fixed = {vkey: volmesh.vertex_coordinates(vkey) for vkey in volmesh.vertices() if TOL.is_close(volmesh.vertex_attribute(vkey, "z"), zmin)}
    assert len(fixed) == 6

    calls = []
    volmesh_planarise(volmesh, kmax=500, fix_vkeys=list(fixed), callback=lambda volmesh, k, args: calls.append(k))

    assert len(calls) < 500
    for vkey in fixed:
        assert volmesh.vertex_coordinates(vkey) == fixed[vkey]
    for hfkey in volmesh.faces():
        for vkey in volmesh.halfface_vertices(hfkey):
            assert distance_point_plane(volmesh.vertex_coordinates(vkey), (volmesh.halfface_center(hfkey), volmesh.halfface_normal(hfkey))) < TOLERANCE_FLAT
    assert volmesh.number_of_cells() == 9
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0


def test_volmesh_planarise_callback_not_callable():
    with pytest.raises(Exception, match="not callable"):
        volmesh_planarise(ForceVolMesh.from_meshgrid(dx=2, nx=2), callback=1)
