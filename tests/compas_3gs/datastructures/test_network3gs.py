import pytest

import compas
from compas.datastructures import Graph
from compas.geometry import Vector
from compas.tolerance import TOL
from compas_3gs.datastructures import Network3gs


@pytest.fixture
def network():
    # a path 0 - 1 - 2 - 3 with edges of length 5, 1 and 3, and a gap in the node keys
    network = Network3gs()
    network.add_node(0, x=0.0, y=0.0, z=0.0)
    network.add_node(1, x=3.0, y=4.0, z=0.0)
    network.add_node(2, x=3.0, y=4.0, z=1.0)
    network.add_node(5, x=3.0, y=1.0, z=1.0)
    network.add_edge(0, 1)
    network.add_edge(1, 2)
    network.add_edge(2, 5)
    return network


def test_constructor():
    network = Network3gs()
    assert isinstance(network, Graph)
    assert network.name == "Network3gs"
    assert network.dual is None
    assert network.attributes == {}
    assert network.default_node_attributes == {"x": 0.0, "y": 0.0, "z": 0.0}
    assert network.default_edge_attributes == {"lmin": 0.0, "lmax": 1e7, "_a": 0.0}
    assert network.number_of_nodes() == 0
    assert Network3gs(name="form").name == "form"


def test_edge_avg_length(network):
    assert [network.edge_length(edge) for edge in network.edges()] == pytest.approx([5.0, 1.0, 3.0])
    assert network.edge_avg_length() == pytest.approx(3.0)


def test_edge_direction_and_edge_vector(network):
    direction = network.edge_direction((0, 1))
    assert isinstance(direction, Vector)
    assert TOL.is_allclose(direction, [0.6, 0.8, 0.0])
    assert TOL.is_allclose(network.edge_direction((1, 0)), [-0.6, -0.8, 0.0])
    assert TOL.is_allclose(network.edge_vector((0, 1)), [3.0, 4.0, 0.0])
    assert TOL.is_allclose(network.edge_vector((2, 5)), [0.0, -3.0, 0.0])
    for edge in network.edges():
        assert TOL.is_close(network.edge_direction(edge).length, 1.0)


def test_vertex_update_xyz_without_fix_attributes(network):
    network.vertex_update_xyz(5, [7.0, 8.0, 9.0])
    assert TOL.is_allclose(network.node_coordinates(5), [7.0, 8.0, 9.0])


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_vertex_update_xyz_with_fix_attributes(network, axis):
    network.node_attribute(1, "xyz"[axis] + "_fix", True)
    network.node_attribute(2, "xyz"[axis] + "_fix", False)

    network.vertex_update_xyz(1, [7.0, 8.0, 9.0])
    expected = [7.0, 8.0, 9.0]
    expected[axis] = [3.0, 4.0, 0.0][axis]
    assert TOL.is_allclose(network.node_coordinates(1), expected)

    network.vertex_update_xyz(2, [7.0, 8.0, 9.0])
    assert TOL.is_allclose(network.node_coordinates(2), [7.0, 8.0, 9.0])

    network.vertex_update_xyz(1, [1.0, 2.0, 3.0], constrained=False)
    assert TOL.is_allclose(network.node_coordinates(1), [1.0, 2.0, 3.0])


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(network, how):
    network.name = "path"
    network.dual = Network3gs()
    network.attributes["note"] = "test"
    network.node_attribute(1, "x_fix", True)
    network.edge_attribute((1, 2), "_a", 12.5)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(network))
    else:
        other = network.copy(copy_guid=True)

    assert type(other) is Network3gs
    assert other.guid == network.guid
    assert other.name == "path"
    assert other.dual is None
    assert other.attributes == {"note": "test"}
    assert other.default_node_attributes == network.default_node_attributes
    assert other.default_edge_attributes == network.default_edge_attributes
    assert list(other.nodes()) == [0, 1, 2, 5]
    assert list(other.edges()) == list(network.edges())
    for node in network.nodes():
        assert TOL.is_allclose(other.node_coordinates(node), network.node_coordinates(node))
        assert other.node_attributes(node) == network.node_attributes(node)
    for edge in network.edges():
        assert other.edge_attributes(edge) == network.edge_attributes(edge)
    assert other.edge_avg_length() == pytest.approx(3.0)
