
from typing import Any

from compas_3gs.datastructures import Network3gs


__all__ = ['CellNetwork']


class CellNetwork(Network3gs):
    """A disjointed, non-polyhedral form and force diagram.

    Each node can hold a cell in its attribute `cell`; the attribute `gfp` holds the global
    force polygon. Both are `None` by default.

    """

    def __init__(self, *args: Any, name: str = "CellNetwork", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)

        self.cells = {}

        self.attributes.update({
            'gfp': None,
        })
        self.default_node_attributes.update({
            'x_fix': False,
            'y_fix': False,
            'z_fix': False,
            'cell': None,
        })
        self.default_edge_attributes.update({
            'target_vector': None,
            'target_length': None,
        })


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
