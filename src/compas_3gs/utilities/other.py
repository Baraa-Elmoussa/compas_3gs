
from collections.abc import Callable
from collections.abc import Hashable
from collections.abc import Sequence

from compas.geometry import dot_vectors

from compas.colors import Color

from compas.datastructures import Graph
from compas.datastructures import VolMesh

from compas_3gs.utilities.topology import pair_uv_to_hf


__all__ = ['get_index_colordict',
           'valuedict_to_colordict',
           'compare_initial_current',

           'get_force_mags',
           'get_force_colors_uv',
           'get_force_colors_hf',

           'print_result']


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   general
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def get_index_colordict(
    key_list: Sequence[Hashable],
    color_scheme: Callable[[float], tuple[int, int, int]] | None = None,
) -> dict[Hashable, tuple[int, int, int]]:
    """Convert a list of keys to a dictionary of key-color pairs based on a color gradient scheme.

    Parameters
    ----------
    key_list
        List of keys.
    color_scheme
        A function that maps a float between 0 and 1 to an RGB tuple (0-255).
        Default is `Color.from_i(value).rgb255` (blue to red).

    Returns
    -------
    dict[Hashable, tuple[int, int, int]]
        Dictionary of key-color pairs, from lowest to highest based on a color gradient scheme.

    """
    c_dict = {}

    if key_list:

        for index, key in enumerate(key_list):

            if len(key_list) == 1:
                value = 1
            else:
                value = float(index) / (len(key_list) - 1)

            color = color_scheme(value) if color_scheme else Color.from_i(max(0.0, min(1.0, value))).rgb255
            c_dict[key] = color

    return c_dict


def valuedict_to_colordict(
    value_dict: dict[Hashable, float],
    color_scheme: Callable[[float], tuple[int, int, int]] | None = None,
) -> dict[Hashable, tuple[int, int, int]]:
    """Convert a dictionary of key-value pairs to a dictionary of key-color pairs.

    Parameters
    ----------
    value_dict
        Dictionary of key-value pairs.
    color_scheme
        A function that maps a float between 0 and 1 to an RGB tuple (0-255).
        Default is `Color.from_i(value).rgb255` (blue to red).

    Returns
    -------
    dict[Hashable, tuple[int, int, int]]
        Dictionary of key-color pairs, from lowest to highest based on a color gradient scheme.

    """
    c_dict = {}

    lb = min(value_dict.values())
    ub = max(value_dict.values())

    for key in value_dict:

        value = (value_dict[key] - lb) / (ub - lb)
        color = color_scheme(value) if color_scheme else Color.from_i(max(0.0, min(1.0, value))).rgb255
        c_dict[key] = color

    return c_dict


def compare_initial_current(
    current_value_dict: dict[Hashable, float],
    initial_value_dict: dict[Hashable, float],
    color_scheme: Callable[[float], tuple[int, int, int]] | None = None,
) -> dict[Hashable, tuple[int, int, int]]:
    """Convert the ratios of current to initial values to a dictionary of key-color pairs.

    Parameters
    ----------
    current_value_dict
        Dictionary of key-value pairs of the current values.
    initial_value_dict
        Dictionary of key-value pairs of the initial values.
    color_scheme
        A function that maps a float to an RGB tuple (0-255).
        Default is `Color.from_i(value).rgb255` (blue to red).

    Returns
    -------
    dict[Hashable, tuple[int, int, int]]
        Dictionary of key-color pairs, for the keys with a current value larger than 0.01.

    """
    color_dict = {}

    for key in current_value_dict:

        current = current_value_dict[key]
        initial = initial_value_dict[key]

        if current > 0.01:

            if initial < 0.01:
                value = 0
            else:
                value = current / initial

            color_dict[key] = color_scheme(value) if color_scheme else Color.from_i(max(0.0, min(1.0, value))).rgb255

    return color_dict


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   force related
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def get_force_mags(volmesh, network):
    """Returns a dictionary of (u,v)-magnitude pairs.
    Negative magnitude means compression, while positive magnitude means tension.
    """

    uv_hf_dict = pair_uv_to_hf(network, volmesh)

    mags = {}

    for u, v in network.edges():
        hfkey = uv_hf_dict[(u, v)]
        edge_vector = network.edge_vector(u, v)
        face_normal = volmesh.halfface_normal(hfkey)
        face_area = volmesh.halfface_area(hfkey)
        dot = dot_vectors(face_normal, edge_vector)

        if dot < 0:
            factor = -1

        if dot > 0:
            factor = 1

        force = face_area * factor * volmesh.attributes['convention']

        mags[(u, v)] = force

    return mags


def get_force_colors_uv(
    volmesh: VolMesh,
    network: Graph,
    gradient: bool = False,
    tol: float = 0.001,
) -> dict[tuple[int, int], tuple[int, int, int]]:
    """Returns a dictionary of (u,v)-color pairs.
    Blue means compression, and red means tension.
    """
    c_dict = {}
    f_dict = get_force_mags(volmesh, network)
    f_range = sorted(f_dict.values())

    c_forces = [x for x in f_range if x < 0]
    t_forces = [x for x in f_range if x > 0]

    for edge in f_dict:
        force = f_dict[edge]
        if force < 0:  # if compression
            color = (0, 0, 255)
            if gradient:
                min_c = abs(c_forces[-1])
                max_c = abs(c_forces[0])
                if min_c != max_c:
                    i = max(0.0, min(1.0, (abs(force) - min_c) / (max_c - min_c)))
                    color = Color(1 - i, 1 - i, 1.0).rgb255

        if force > 0:  # if tension
            color = (255, 0, 0)
            if gradient:
                min_t = t_forces[0]
                max_t = t_forces[-1]
                if min_t != max_t:
                    i = max(0.0, min(1.0, (force - min_t) / (max_t - min_t)))
                    color = Color(1.0, 1 - i, 1 - i).rgb255

        if abs(force) < tol:  # if close to zero
            color = (255, 255, 255)

        c_dict[edge] = color

    return c_dict


def get_force_colors_hf(
    volmesh: VolMesh,
    network: Graph,
    uv_c_dict: dict[tuple[int, int], tuple[int, int, int]] | None = None,
    gradient: bool = False,
    boundary: bool = False,
    tol: float = 0.001,
) -> dict[int, tuple[int, int, int]]:
    """Returns a dictionary of hfkey-color pairs.
    """
    uv_hf_dict = pair_uv_to_hf(network, volmesh)

    uv_c_dict = uv_c_dict or get_force_colors_uv(volmesh,
                                                 network,
                                                 gradient=gradient)

    # interior halffaces -------------------------------------------------------
    hf_c_dict = {}
    for uv in uv_c_dict:
        u_hfkey = uv_hf_dict[uv]
        v_hfkey = volmesh.halfface_opposite_halfface(u_hfkey)

        hf_c_dict[u_hfkey] = uv_c_dict[uv]
        hf_c_dict[v_hfkey] = uv_c_dict[uv]

    # boundary halffaces -------------------------------------------------------
    if boundary:
        b_hfkeys = volmesh.halffaces_on_boundaries()
        b_hf_areas = {hfkey: volmesh.halfface_area(hfkey) for hfkey in b_hfkeys}
        b_hf_c_dict = valuedict_to_colordict(b_hf_areas,
                                             color_scheme=lambda i: Color(1 - max(0.0, min(1.0, i)), 1.0, 1 - max(0.0, min(1.0, i))).rgb255)
        for hfkey in b_hfkeys:
            hf_color = (0, 255, 0)

            if gradient:
                hf_color = b_hf_c_dict[hfkey]
            hf_c_dict[hfkey] = hf_color

    return hf_c_dict


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   other
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def print_result(name, k, deviation):
    """Prints the result of an algorithm.

    Parameters
    ----------
    name : string
        name of the algorithm
    k : int
        number of iterations
    deviation : float
        deviation

    """
    name = str(name)

    print('===================================================================')
    print('')
    print(name, 'ended after', k, 'iterations.')
    print('')
    print('Max deviation :', deviation)
    print('')
    print('===================================================================')


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
