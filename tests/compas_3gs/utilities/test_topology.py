import json
import os

import pytest

from compas.geometry import angle_vectors
from compas.geometry import dot_vectors
from compas.tolerance import TOL
from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.datastructures import VolMesh3gs
from compas_3gs.utilities import pair_hf_to_uv
from compas_3gs.utilities import pair_uv_to_hf


@pytest.fixture(params=["meshgrid", "honeycomb"])
def diagrams(request):
    if request.param == "meshgrid":
        volmesh = VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)
    else:
        with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
            data = json.load(f)
        volmesh = VolMesh3gs.from_vertices_and_cells(data["vertices"], data["cells"])
    network = volmesh_dual_network(volmesh)
    return volmesh, network


def test_pair_uv_to_hf(diagrams):
    volmesh, network = diagrams
    uv_hf = pair_uv_to_hf(network, volmesh)

    assert sorted(uv_hf) == sorted(network.edges())
    for (u, v), hf in uv_hf.items():
        assert volmesh.halfface_cell(hf) == u
        assert volmesh.halfface_opposite_cell(hf) == v
        # the halfface of cell u points towards cell v
        assert dot_vectors(volmesh.halfface_normal(hf), network.edge_vector((u, v))) > 0


def test_pair_hf_to_uv(diagrams):
    volmesh, network = diagrams
    hf_uv = pair_hf_to_uv(volmesh, network)

    # every interior halfface, once
    interior = [hf for hf in volmesh.halffaces() if not volmesh.is_halfface_on_boundary(hf)]
    assert len(hf_uv) == 2 * network.number_of_edges()
    assert sorted(hf_uv) == sorted(interior)
    for hf, (u, v) in hf_uv.items():
        assert volmesh.halfface_cell(hf) == u
        assert volmesh.halfface_opposite_cell(hf) == v


def test_pair_hf_to_uv_inverse(diagrams):
    volmesh, network = diagrams
    uv_hf = pair_uv_to_hf(network, volmesh)
    hf_uv = pair_hf_to_uv(volmesh, network)

    for uv, hf in uv_hf.items():
        assert hf_uv[hf] == uv
        assert hf_uv[volmesh.halfface_opposite_halfface(hf)] == uv[::-1]
    for hf, (u, v) in hf_uv.items():
        if network.has_edge((u, v)):
            assert uv_hf[(u, v)] == hf
        else:
            assert uv_hf[(v, u)] == volmesh.halfface_opposite_halfface(hf)


def test_pair_uv_to_hf_meshgrid_normals():
    volmesh = VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)
    network = volmesh_dual_network(volmesh)
    uv_hf = pair_uv_to_hf(network, volmesh)

    # 8 cubes, 12 shared faces; every dual edge is parallel to the normal of its halfface
    assert len(uv_hf) == 12
    for uv, hf in uv_hf.items():
        assert angle_vectors(volmesh.halfface_normal(hf), network.edge_vector(uv)) < TOL.angular
        assert TOL.is_close(volmesh.halfface_area(hf), 0.25)


def test_pair_uv_to_hf_not_neighbours():
    volmesh = VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)
    network = volmesh_dual_network(volmesh)
    # cells 0 and 7 are opposite corners of the grid
    assert not network.has_edge((0, 7), directed=False)
    network.add_edge(0, 7)

    with pytest.raises(KeyError):
        pair_uv_to_hf(network, volmesh)
    with pytest.raises(KeyError):
        pair_hf_to_uv(volmesh, network)
