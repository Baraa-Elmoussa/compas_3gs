import pytest

from compas.datastructures import Graph
from compas.datastructures import Mesh
from compas.datastructures import VolMesh
from compas.geometry import Box
from compas.geometry import centroid_points
from compas.geometry import volume_polyhedron
from compas_3gs.datastructures import VolMesh3gs
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.operations import point_reflection
from compas_3gs.operations import volmesh_cell_split_vertex
from compas_3gs.operations import volmesh_cell_subdivide_barycentric
from compas_3gs.operations import volmesh_halfface_extrude
from compas_3gs.operations import volmesh_halfface_pinch
from compas_3gs.operations import volmesh_halfface_subdivide
from compas_3gs.operations import volmesh_merge_adjacent_halffaces
from compas_3gs.operations import volmesh_vertex_lift
from compas_3gs.operations import volmesh_vertex_merge

# `VolMesh.is_valid()` is not implemented in compas 2.15.1. A volmesh is checked with: every cell closed,
# every halfface pointing out of its cell (signed volume > 0; `Mesh.volume` returns the absolute value),
# and every halfface that shares its vertices with another one paired with it as opposites.


@pytest.fixture
def meshgrid():
    # 2 x 2 x 2 cubes of 0.5, volume 1
    return VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)


@pytest.fixture
def lshape():
    # 2 x 2 x 1 unit cubes without the cell at x, y in [1, 2]: a reflex edge at x = y = 1, volume 3
    volmesh = VolMesh3gs.from_meshgrid(2, 2, 1, 2, 2, 1)
    volmesh.delete_cell(3)
    volmesh.remove_unused_vertices()
    return volmesh


@pytest.fixture(params=["meshgrid", "lshape"])
def volmesh(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def slab():
    # 3 x 3 x 1 unit cubes, volume 9
    return VolMesh3gs.from_meshgrid(3, 3, 1, 3, 3, 1)


# ==============================================================================
# global operations
# ==============================================================================


def test_point_reflection_graph():
    graph = Graph()
    for x, y, z in [[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 4]]:
        graph.add_node(x=x, y=y, z=z)
    graph.add_edge(0, 1)
    graph.add_edge(0, 3)
    center = [0.5, 0.25, 1.0]
    xyz = {node: graph.node_coordinates(node) for node in graph.nodes()}
    vectors = {edge: graph.edge_vector(edge) for edge in graph.edges()}

    assert point_reflection(graph) is None

    for node in graph.nodes():
        assert graph.node_coordinates(node) == pytest.approx([2 * c - p for c, p in zip(center, xyz[node])])
    for edge in graph.edges():
        assert list(graph.edge_vector(edge)) == pytest.approx([-c for c in vectors[edge]])


def test_point_reflection_mesh():
    mesh = Mesh.from_shape(Box(1))
    mesh.vertex_attributes(6, "xyz", [0.9, 0.7, 1.1])
    center = centroid_points(mesh.vertices_attributes("xyz"))
    xyz = {vertex: mesh.vertex_coordinates(vertex) for vertex in mesh.vertices()}
    normals = {face: mesh.face_normal(face) for face in mesh.faces()}
    volume = volume_polyhedron(mesh.to_vertices_and_faces())
    assert volume > 0

    point_reflection(mesh)

    for vertex in mesh.vertices():
        assert mesh.vertex_coordinates(vertex) == pytest.approx([2 * c - p for c, p in zip(center, xyz[vertex])])
    # every face normal keeps its direction, so the closed mesh turns inside out
    for face in mesh.faces():
        assert list(mesh.face_normal(face)) == pytest.approx(list(normals[face]))
    assert mesh.is_closed()
    assert volume_polyhedron(mesh.to_vertices_and_faces()) == pytest.approx(-volume)


@pytest.mark.parametrize("cls", [VolMesh, ForceVolMesh])
def test_point_reflection_volmesh(cls):
    volmesh = cls.from_meshgrid(1, 1, 1, 2, 2, 2)
    volmesh.vertex_attributes(26, "xyz", [1.3, 1.1, 1.2])
    center = centroid_points(volmesh.vertices_attributes("xyz"))
    xyz = {vertex: volmesh.vertex_coordinates(vertex) for vertex in volmesh.vertices()}
    normals = {hf: volmesh.halfface_normal(hf) for hf in volmesh.halffaces()}
    volumes = {cell: volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()}
    assert all(volume > 0 for volume in volumes.values())

    assert point_reflection(volmesh) is None

    for vertex in volmesh.vertices():
        assert volmesh.vertex_coordinates(vertex) == pytest.approx([2 * c - p for c, p in zip(center, xyz[vertex])])
    assert centroid_points(volmesh.vertices_attributes("xyz")) == pytest.approx(center)
    for hf in volmesh.halffaces():
        assert list(volmesh.halfface_normal(hf)) == pytest.approx(list(normals[hf]))
    # the cells turn inside out on purpose (compression / tension switch); the topology is unchanged
    assert volmesh.number_of_cells() == 8
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) == pytest.approx(-volumes[cell])
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


# ==============================================================================
# vertex operations
# ==============================================================================


def test_volmesh_vertex_lift(meshgrid):
    volmesh = meshgrid
    vertex = [v for v in volmesh.vertices() if volmesh.vertex_coordinates(v) == [0.5, 0.5, 1.0]][0]
    halffaces = [hf for hf in volmesh.vertex_halffaces(vertex) if volmesh.is_halfface_on_boundary(hf)]
    assert len(halffaces) == 4

    result = volmesh_vertex_lift(volmesh, vertex, [0.5, 0.5, 1.5], halffaces)

    assert result is volmesh
    # one pyramid per halfface, together a pyramid of height 0.5 on the unit top face
    assert volmesh.number_of_cells() == 12
    assert volmesh.number_of_vertices() == 28
    for hf in halffaces:
        assert not volmesh.is_halfface_on_boundary(hf)
        assert len(volmesh.cell_faces(volmesh.halfface_opposite_cell(hf))) == 5
    assert sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()) == pytest.approx(1 + 1 / 6)
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) > 0
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


def test_volmesh_vertex_lift_interior(meshgrid):
    vertex = [v for v in meshgrid.vertices() if meshgrid.vertex_coordinates(v) == [0.5, 0.5, 0.5]][0]
    with pytest.raises(Exception, match="interior"):
        volmesh_vertex_lift(meshgrid, vertex, [0.5, 0.5, 1.5], [])


def test_volmesh_vertex_merge(slab):
    volmesh = slab
    xyz = {tuple(volmesh.vertex_coordinates(v)): v for v in volmesh.vertices()}
    u, v = xyz[(1.0, 1.0, 1.0)], xyz[(1.0, 2.0, 1.0)]
    # u and v are the first and last vertex of a halfface: the merged face wraps around
    assert any([volmesh.halfface_vertices(hf)[0], volmesh.halfface_vertices(hf)[-1]] == [u, v] for hf in volmesh.halffaces())

    result = volmesh_vertex_merge(volmesh, [u, v], [1.0, 1.5, 1.0])

    assert result is volmesh
    assert volmesh.number_of_cells() == 9
    assert volmesh.number_of_vertices() == 31
    assert not volmesh.has_vertex(u) and not volmesh.has_vertex(v)
    w = [vertex for vertex in volmesh.vertices() if volmesh.vertex_coordinates(vertex) == [1.0, 1.5, 1.0]]
    assert len(w) == 1
    # the two top faces and the vertical face on the edge u-v become triangles
    for hf in volmesh.halffaces():
        vertices = volmesh.halfface_vertices(hf)
        assert len(set(vertices)) == len(vertices)
    assert len([face for face in volmesh.faces() if len(volmesh.halfface_vertices(face)) == 3]) == 3
    # the top stays in the plane z = 1, so the volume does not change
    assert sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()) == pytest.approx(9)
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) > 0
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


def test_volmesh_vertex_merge_attributes(slab):
    volmesh = slab
    xyz = {tuple(volmesh.vertex_coordinates(v)): v for v in volmesh.vertices()}
    merged = {xyz[(1.0, 1.0, 1.0)], xyz[(1.0, 2.0, 1.0)]}
    for cell in volmesh.cells():
        volmesh.cell_attribute(cell, "tag", cell)
    faces = {}
    for face in volmesh.faces():
        vertices = volmesh.halfface_vertices(face)
        volmesh.face_attribute(face, "tag", str(sorted(vertices)))
        faces[frozenset("w" if key in merged else key for key in vertices)] = str(sorted(vertices))
    edges = {}
    for edge in volmesh.edges():
        volmesh.edge_attribute(edge, "tag", str(sorted(edge)))
        if set(edge) != merged:
            edges[frozenset("w" if key in merged else key for key in edge)] = str(sorted(edge))

    volmesh_vertex_merge(volmesh, list(merged), [1.0, 1.5, 1.0])

    w = [vertex for vertex in volmesh.vertices() if volmesh.vertex_coordinates(vertex) == [1.0, 1.5, 1.0]][0]
    assert {cell: volmesh.cell_attribute(cell, "tag") for cell in volmesh.cells()} == {cell: cell for cell in range(9)}
    for hf in volmesh.halffaces():
        assert volmesh.face_attribute(hf, "tag") == faces[frozenset("w" if key == w else key for key in volmesh.halfface_vertices(hf))]
    assert len(list(volmesh.edges())) == len(edges)
    for edge in volmesh.edges():
        assert volmesh.edge_attribute(edge, "tag") == edges[frozenset("w" if key == w else key for key in edge)]


@pytest.mark.parametrize("points", [[(1.0, 1.0, 1.0), (2.0, 2.0, 1.0)], [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (2.0, 2.0, 1.0)]])
def test_volmesh_vertex_merge_illegal(slab, points):
    # opposite corners of a quad, or three of its four corners
    xyz = {tuple(slab.vertex_coordinates(v)): v for v in slab.vertices()}
    with pytest.raises(Exception, match="Illegal merge"):
        volmesh_vertex_merge(slab, [xyz[point] for point in points], [1.5, 1.5, 1.0])


# ==============================================================================
# halfface operations
# ==============================================================================


def test_volmesh_halfface_pinch(volmesh):
    volume = sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells())
    cells = volmesh.number_of_cells()
    hf = [hf for hf in sorted(volmesh.halffaces_on_boundaries()) if volmesh.halfface_normal(hf)[2] > 0.9][0]
    x, y, z = volmesh.halfface_center(hf)
    area = volmesh.halfface_area(hf)

    result = volmesh_halfface_pinch(volmesh, hf, [x, y, z + 0.3])

    assert result is volmesh
    assert volmesh.number_of_cells() == cells + 1
    assert not volmesh.is_halfface_on_boundary(hf)
    pyramid = volmesh.halfface_opposite_cell(hf)
    assert len(volmesh.cell_faces(pyramid)) == 5
    assert [x, y, z + 0.3] in [volmesh.vertex_coordinates(vertex) for vertex in volmesh.cell_vertices(pyramid)]
    assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(pyramid)) == pytest.approx(area * 0.3 / 3)
    assert sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()) == pytest.approx(volume + area * 0.3 / 3)
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) > 0
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


def test_volmesh_merge_adjacent_halffaces_reflex_edge(lshape):
    volmesh = lshape
    # the two boundary halffaces that face the missing cell meet at the reflex edge
    notch = [hf for hf in volmesh.halffaces_on_boundaries() if list(volmesh.halfface_center(hf)) in ([1.0, 1.5, 0.5], [1.5, 1.0, 0.5])]
    assert len(notch) == 2

    result = volmesh_merge_adjacent_halffaces(volmesh, notch)

    assert result is volmesh
    assert volmesh.number_of_cells() == 4
    # the new cell is the convex hull of the 6 vertices: half of the missing cube, a triangular prism
    cell = volmesh.halfface_opposite_cell(notch[0])
    assert cell == volmesh.halfface_opposite_cell(notch[1])
    assert cell not in (0, 1, 2)
    assert sorted(len(volmesh.halfface_vertices(hf)) for hf in volmesh.cell_faces(cell)) == [3, 3, 4, 4, 4]
    assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) == pytest.approx(0.5)
    assert sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()) == pytest.approx(3.5)
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) > 0
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


def test_volmesh_merge_adjacent_halffaces_errors(lshape):
    interior = [hf for hf in lshape.halffaces() if not lshape.is_halfface_on_boundary(hf)][0]
    with pytest.raises(ValueError, match="interior"):
        volmesh_merge_adjacent_halffaces(lshape, [interior])

    # the bottom of cell 0 and the top of cell 2 share no edge
    bottom = [hf for hf in lshape.cell_faces(0) if lshape.halfface_normal(hf)[2] < -0.9][0]
    top = [hf for hf in lshape.cell_faces(2) if lshape.halfface_normal(hf)[2] > 0.9][0]
    with pytest.raises(ValueError, match="not chained"):
        volmesh_merge_adjacent_halffaces(lshape, [bottom, top])
    assert lshape.number_of_cells() == 3


# ==============================================================================
# cell operations
# ==============================================================================


def test_volmesh_cell_subdivide_barycentric(volmesh):
    volume = sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells())
    cells = volmesh.number_of_cells()
    cell = next(iter(volmesh.cells()))
    center = list(volmesh.cell_center(cell))
    cell_volume = volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell))

    new_cells = volmesh_cell_subdivide_barycentric(volmesh, cell)

    # one pyramid per face of the cube, each a sixth of its volume, with the apex at the cube centre
    assert len(new_cells) == 6
    assert cell not in volmesh.cells()
    assert volmesh.number_of_cells() == cells + 5
    for new_cell in new_cells:
        assert len(volmesh.cell_faces(new_cell)) == 5
        assert center in [volmesh.vertex_coordinates(vertex) for vertex in volmesh.cell_vertices(new_cell)]
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(new_cell)) == pytest.approx(cell_volume / 6)
    assert sum(volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) for cell in volmesh.cells()) == pytest.approx(volume)
    for cell in volmesh.cells():
        assert volmesh.cell_to_mesh(cell).is_closed()
        assert volume_polyhedron(volmesh.cell_to_vertices_and_faces(cell)) > 0
    for hf in volmesh.halffaces():
        twins = [other for other in volmesh.halffaces() if other != hf and set(volmesh.halfface_vertices(other)) == set(volmesh.halfface_vertices(hf))]
        assert twins == ([] if volmesh.is_halfface_on_boundary(hf) else [volmesh.halfface_opposite_halfface(hf)])


def test_volmesh_cell_subdivide_barycentric_attributes(meshgrid):
    volmesh = meshgrid
    volmesh.cell_attribute(0, "tag", "old")
    for hf in volmesh.cell_faces(0):
        volmesh.face_attribute(hf, "tag", str(sorted(volmesh.halfface_vertices(hf))))
    for edge in volmesh.cell_edges(0):
        volmesh.edge_attribute(edge, "tag", str(sorted(edge)))
    edges = {frozenset(edge) for edge in volmesh.cell_edges(0)}

    new_cells = volmesh_cell_subdivide_barycentric(volmesh, 0)

    for cell in new_cells:
        assert volmesh.cell_attribute(cell, "tag") is None
        base = [hf for hf in volmesh.cell_faces(cell) if len(volmesh.halfface_vertices(hf)) == 4][0]
        assert volmesh.face_attribute(base, "tag") == str(sorted(volmesh.halfface_vertices(base)))
        for edge in volmesh.cell_edges(cell):
            assert volmesh.edge_attribute(edge, "tag") == (str(sorted(edge)) if frozenset(edge) in edges else None)


@pytest.mark.parametrize(
    "function, args",
    [
        (volmesh_cell_split_vertex, (0, 0, 1)),
        (volmesh_halfface_extrude, (0,)),
        (volmesh_halfface_subdivide, (0,)),
    ],
)
def test_volmesh_not_implemented(meshgrid, function, args):
    with pytest.raises(NotImplementedError):
        function(meshgrid, *args)
