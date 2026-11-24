
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from compas.datastructures import Graph
    from compas.datastructures import VolMesh


__all__ = ['pair_hf_to_uv',
           'pair_uv_to_hf']


def pair_hf_to_uv(volmesh: "VolMesh", network: "Graph") -> dict[int, tuple[int, int]]:
    """Pairs the directed halffaces of a volmesh to the corresponding edge (u, v) of the dual network.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    network
        A network object representing a polyhedral form diagram.

    Returns
    -------
    dict[int, tuple[int, int]]
        A dictionary of halfface-edge pairs: the halfface of cell `u` gives `(u, v)`,
        its opposite halfface, of cell `v`, gives `(v, u)`.

    Notes
    -----
    `u_hfkey` is an interior halfface of the volmesh that belongs to cell `u` and points to cell `v`.
    In other words, its pair (or opposite) halfface belongs to cell `v`.

    """
    hf_uv_dict = {}

    for u, v in network.edges():
        u_hfkey, v_hfkey = volmesh.cell_pair_halffaces(u, v)
        hf_uv_dict[u_hfkey] = (u, v)
        hf_uv_dict[v_hfkey] = (v, u)

    return hf_uv_dict


def pair_uv_to_hf(network: "Graph", volmesh: "VolMesh") -> dict[tuple[int, int], int]:
    """Pairs the directed edges (u, v) of a network to the corresponding halffaces of the dual volmesh.

    Parameters
    ----------
    network
        A network object representing a polyhedral form diagram.
    volmesh
        A volmesh object representing a polyhedral force diagram.

    Returns
    -------
    dict[tuple[int, int], int]
        A dictionary of edge-halfface pairs: `(u, v)` gives the halfface of cell `u` that points to cell `v`.

    """
    uv_hf_dict = {}

    for u, v in network.edges():
        u_hfkey, v_hfkey = volmesh.cell_pair_halffaces(u, v)
        uv_hf_dict[(u, v)] = u_hfkey

    return uv_hf_dict


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == '__main__':
    pass
