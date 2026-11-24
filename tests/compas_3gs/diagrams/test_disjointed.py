import pytest

import compas
from compas.geometry import Box
from compas.geometry import Frame
from compas.tolerance import TOL
from compas_3gs.datastructures import Network3gs
from compas_3gs.diagrams import Cell
from compas_3gs.diagrams import CellNetwork


@pytest.fixture
def network():
    network = CellNetwork()
    network.add_node(0, x=0.0, y=0.0, z=0.0)
    network.add_node(1, x=2.0, y=0.0, z=0.0)
    network.add_node(3, x=2.0, y=2.0, z=0.0)
    network.add_edge(0, 1)
    network.add_edge(1, 3)
    return network


def test_constructor():
    network = CellNetwork()
    assert isinstance(network, Network3gs)
    assert network.name == "CellNetwork"
    assert network.dual is None
    assert network.cells == {}
    assert network.attributes == {"gfp": None}
    assert network.default_node_attributes == {"x": 0.0, "y": 0.0, "z": 0.0, "x_fix": False, "y_fix": False, "z_fix": False, "cell": None}
    for name, value in {"target_vector": None, "target_length": None}.items():
        assert network.default_edge_attributes[name] == value
    assert CellNetwork(name="cells").name == "cells"


def test_nodes_get_none_and_do_not_share_a_cell(network):
    for node in network.nodes():
        assert network.node_attribute(node, "cell") is None

    cell = Cell.from_shape(Box(1))
    network.node_attribute(0, "cell", cell)
    assert network.node_attribute(0, "cell") is cell
    assert network.node_attribute(1, "cell") is None
    assert network.node_attribute(3, "cell") is None
    assert network.node_attribute(network.add_node(), "cell") is None

    network.node_attribute(1, "cell", Cell.from_shape(Box(2)))
    assert network.node_attribute(1, "cell") is not network.node_attribute(0, "cell")
    assert network.node_attribute(0, "cell").number_of_vertices() == 8


def test_instances_do_not_share_state():
    a = CellNetwork()
    b = CellNetwork()
    a.cells[0] = Cell()
    a.attributes["gfp"] = Cell()
    assert b.cells == {}
    assert b.attributes["gfp"] is None
    assert CellNetwork().default_node_attributes["cell"] is None


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(network, how):
    network.name = "disjointed"
    network.attributes["gfp"] = Cell.from_shape(Box(1))
    network.node_attribute(0, "cell", Cell.from_shape(Box(0.5, frame=Frame([0.0, 0.0, 0.0]))))
    network.node_attribute(1, "cell", Cell.from_shape(Box(0.5, frame=Frame([2.0, 0.0, 0.0]))))
    network.node_attribute(1, "x_fix", True)
    network.edge_attribute((0, 1), "target_length", 2.0)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(network))
    else:
        other = network.copy(copy_guid=True)

    assert type(other) is CellNetwork
    assert other.guid == network.guid
    assert other.name == "disjointed"
    assert other.dual is None
    assert other.attributes.keys() == {"gfp"}
    assert other.default_node_attributes == network.default_node_attributes
    assert other.default_edge_attributes == network.default_edge_attributes
    assert list(other.nodes()) == [0, 1, 3]
    assert list(other.edges()) == list(network.edges())
    for node in network.nodes():
        assert TOL.is_allclose(other.node_coordinates(node), network.node_coordinates(node))
    assert other.node_attribute(1, "x_fix") is True
    assert other.node_attribute(3, "cell") is None
    for edge in network.edges():
        assert other.edge_attributes(edge) == network.edge_attributes(edge)

    cells = [(other.attributes["gfp"], network.attributes["gfp"])]
    cells += [(other.node_attribute(node, "cell"), network.node_attribute(node, "cell")) for node in (0, 1)]
    for copied, cell in cells:
        assert type(copied) is Cell
        assert copied is not cell
        assert list(copied.faces()) == list(cell.faces())
        for face in cell.faces():
            assert copied.face_vertices(face) == cell.face_vertices(face)
        for vertex in cell.vertices():
            assert TOL.is_allclose(copied.vertex_coordinates(vertex), cell.vertex_coordinates(vertex))
        assert TOL.is_allclose(copied.cell_center(), cell.cell_center())
    assert other.node_attribute(0, "cell") is not other.node_attribute(1, "cell")
