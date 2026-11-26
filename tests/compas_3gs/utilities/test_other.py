import pytest

from compas.tolerance import TOL
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import compare_initial_current
from compas_3gs.utilities import get_force_colors_hf
from compas_3gs.utilities import get_force_colors_uv
from compas_3gs.utilities import get_force_mags
from compas_3gs.utilities import get_index_colordict
from compas_3gs.utilities import print_result
from compas_3gs.utilities import valuedict_to_colordict

BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
RED = (255, 0, 0)
WHITE = (255, 255, 255)


@pytest.fixture
def diagrams():
    # three unit cubes in a row along x; the face between cells 1 and 2 (x = 2) is enlarged to 2 x 2,
    # so the two interior faces have areas 1 and 4, and the boundary faces areas 1 and sqrt(1.25) * 1.5
    volmesh = ForceVolMesh.from_meshgrid(3, 1, 1, 3, 1, 1)
    for vertex in volmesh.vertices():
        x, y, z = volmesh.vertex_coordinates(vertex)
        if x == 2:
            volmesh.vertex_attributes(vertex, "yz", [2 * y - 0.5, 2 * z - 0.5])
    network = volmesh_dual_network(volmesh, FormNetwork)
    return volmesh, network


# ==============================================================================
# general
# ==============================================================================


def test_get_index_colordict():
    keys = ["a", "b", "c", "d", "e"]
    colors = get_index_colordict(keys)

    assert list(colors) == keys
    for color in colors.values():
        assert isinstance(color, tuple) and len(color) == 3
        assert all(isinstance(c, int) and 0 <= c <= 255 for c in color)
    assert colors["a"] == BLUE
    assert colors["c"] == GREEN
    assert colors["e"] == RED


def test_get_index_colordict_single_key():
    assert get_index_colordict(["a"]) == {"a": RED}


def test_get_index_colordict_empty():
    assert get_index_colordict([]) == {}


def test_get_index_colordict_color_scheme():
    colors = get_index_colordict([10, 20, 30], color_scheme=lambda i: (round(255 * i), 0, 0))
    assert colors == {10: (0, 0, 0), 20: (128, 0, 0), 30: (255, 0, 0)}


def test_valuedict_to_colordict():
    colors = valuedict_to_colordict({"a": -2.0, "b": 0.0, "c": 2.0, "d": 1.0})

    for color in colors.values():
        assert isinstance(color, tuple) and len(color) == 3
        assert all(isinstance(c, int) and 0 <= c <= 255 for c in color)
    assert colors["a"] == BLUE
    assert colors["b"] == GREEN
    assert colors["c"] == RED
    assert colors["d"] == (255, 255, 0)


def test_valuedict_to_colordict_color_scheme():
    colors = valuedict_to_colordict({"a": 3.0, "b": 4.0, "c": 7.0}, color_scheme=lambda i: (round(255 * i), 0, 0))
    assert colors == {"a": (0, 0, 0), "b": (64, 0, 0), "c": (255, 0, 0)}


def test_valuedict_to_colordict_equal_values():
    colors = valuedict_to_colordict({"a": 1.0, "b": 1.0, "c": 1.0})

    assert list(colors) == ["a", "b", "c"]
    # all of them the colour of the lowest value
    assert set(colors.values()) == {valuedict_to_colordict({"low": 0.0, "high": 1.0})["low"]}


def test_compare_initial_current():
    current = {"zero_initial": 1.0, "half": 0.5, "same": 2.0, "double": 2.0, "tiny": 0.005}
    initial = {"zero_initial": 0.0, "half": 1.0, "same": 2.0, "double": 1.0, "tiny": 1.0}
    colors = compare_initial_current(current, initial)

    # keys with a current value of at most 0.01 are left out
    assert set(colors) == {"zero_initial", "half", "same", "double"}
    assert colors["zero_initial"] == BLUE
    assert colors["half"] == GREEN
    assert colors["same"] == RED
    # a ratio above 1 is clamped to 1
    assert colors["double"] == RED


# ==============================================================================
# force related
# ==============================================================================


def test_get_force_mags(diagrams):
    volmesh, network = diagrams

    # every form edge (u, v) runs along the outward normal of its halfface of cell u: compression
    assert get_force_mags(volmesh, network) == {(0, 1): pytest.approx(-1.0), (1, 2): pytest.approx(-4.0)}

    volmesh.attributes["convention"] = 1
    assert get_force_mags(volmesh, network) == {(0, 1): pytest.approx(1.0), (1, 2): pytest.approx(4.0)}


def test_get_force_mags_reversed_edge(diagrams):
    volmesh, network = diagrams

    # node 2 behind node 1: the edge (1, 2) runs against the normal, tension
    network.node_attribute(2, "x", 1.0)
    assert get_force_mags(volmesh, network) == {(0, 1): pytest.approx(-1.0), (1, 2): pytest.approx(4.0)}


def test_get_force_mags_edge_in_the_plane_of_its_face(diagrams):
    volmesh, network = diagrams

    # node 2 beside node 1: the edge (1, 2) is perpendicular to the normal of its face (x = 2)
    x, y, z = network.node_coordinates(1)
    network.node_attributes(2, "xyz", [x, y + 1.0, z])
    assert get_force_mags(volmesh, network) == {(0, 1): pytest.approx(-1.0), (1, 2): 0.0}


def test_get_force_colors_uv(diagrams):
    volmesh, network = diagrams

    assert get_force_colors_uv(volmesh, network) == {(0, 1): BLUE, (1, 2): BLUE}

    network.node_attribute(2, "x", 1.0)
    assert get_force_colors_uv(volmesh, network) == {(0, 1): BLUE, (1, 2): RED}


def test_get_force_colors_uv_tol(diagrams):
    volmesh, network = diagrams

    assert get_force_colors_uv(volmesh, network, tol=2.0) == {(0, 1): WHITE, (1, 2): BLUE}


def test_get_force_colors_uv_gradient(diagrams):
    volmesh, network = diagrams

    # white for the smallest magnitude, the full colour for the largest
    assert get_force_colors_uv(volmesh, network, gradient=True) == {(0, 1): WHITE, (1, 2): BLUE}

    volmesh.attributes["convention"] = 1
    assert get_force_colors_uv(volmesh, network, gradient=True) == {(0, 1): WHITE, (1, 2): RED}


def test_get_force_colors_hf(diagrams):
    volmesh, network = diagrams
    colors = get_force_colors_hf(volmesh, network)

    interior = [hf for hf in volmesh.halffaces() if not volmesh.is_halfface_on_boundary(hf)]
    assert sorted(colors) == sorted(interior)
    for u, v in network.edges():
        hf, opposite = volmesh.cell_pair_halffaces(u, v)
        assert colors[hf] == BLUE
        assert colors[opposite] == BLUE


def test_get_force_colors_hf_uv_colors(diagrams):
    volmesh, network = diagrams
    colors = get_force_colors_hf(volmesh, network, uv_c_dict={(0, 1): GREEN, (1, 2): RED})

    for (u, v), color in [((0, 1), GREEN), ((1, 2), RED)]:
        hf, opposite = volmesh.cell_pair_halffaces(u, v)
        assert colors[hf] == color
        assert colors[opposite] == color


def test_get_force_colors_hf_boundary(diagrams):
    volmesh, network = diagrams
    boundary = volmesh.halffaces_on_boundaries()

    colors = get_force_colors_hf(volmesh, network, boundary=True)
    assert sorted(colors) == sorted(volmesh.halffaces())
    assert all(colors[hf] == GREEN for hf in boundary)

    # white for the smallest boundary area, green for the largest
    colors = get_force_colors_hf(volmesh, network, gradient=True, boundary=True)
    areas = {hf: volmesh.halfface_area(hf) for hf in boundary}
    smallest = min(areas.values())
    largest = max(areas.values())
    assert smallest == pytest.approx(1.0)
    assert largest == pytest.approx(1.5 * 1.25**0.5)
    for hf, area in areas.items():
        assert colors[hf] == (WHITE if TOL.is_close(area, smallest) else GREEN)


# ==============================================================================
# other
# ==============================================================================


def test_print_result(capsys):
    print_result("Planarisation", 12, 0.001)
    out = capsys.readouterr().out
    assert "Planarisation ended after 12 iterations." in out
    assert "Max deviation : 0.001" in out
