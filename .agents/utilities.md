# Utilities

`from compas_3gs.utilities import ...`

## Forces and colours of a form–force pair

| Function | Returns |
|---|---|
| `get_force_mags(volmesh, network)` | `{(u, v): force}` from the face areas, with the sign from `volmesh.attributes["convention"]` (negative: compression); 0 for a member in the plane of its face |
| `get_force_colors_uv(volmesh, network, gradient=False, tol=0.001)` | `{(u, v): (r, g, b)}`, 0–255: blue compression, red tension, white zero (the 3GS convention, which `GSViewer` follows) |
| `get_force_colors_hf(volmesh, network, uv_c_dict=None, gradient=False, boundary=False, tol=0.001)` | `{halfface: (r, g, b)}` |
| `pair_hf_to_uv(volmesh, network)`, `pair_uv_to_hf(network, volmesh)` | halfface ↔ edge of the dual |
| `get_index_colordict(key_list, color_scheme=None)` | a colour per key along a gradient |
| `valuedict_to_colordict(value_dict, color_scheme=None)` | a colour per value, from the lowest to the highest; all the lowest colour if all values are equal |
| `compare_initial_current(current_value_dict, initial_value_dict, color_scheme=None)` | colour by the ratio of current to initial values |

## Geometry

| Function | Returns |
|---|---|
| `cell_face_flatness(cell)`, `volmesh_face_flatness(volmesh)` | `{face: flatness}` |
| `cell_face_areaness(cell, target_areas)`, `volmesh_face_areaness(volmesh, target_areas)` | `{face: abs(target area - area)}` for the faces with a target |
| `polygon_normal_oriented(polygon, unitized=True)`, `polygon_area_oriented(polygon)`, `polygon_area_footprint(polygon)` | normal and areas of a polygon |
| `scale_polygon(polygon, scale)` | the polygon scaled about its centroid |
| `resultant_vector(vectors, locations)` | `(resultant, location)` of vectors acting at locations |
| `datastructure_centroid(datastructure)` | vertex centroid of a mesh or volmesh; for a graph use `centroid_points(graph.nodes_attributes("xyz"))` |
| `print_result(name, k, deviation)` | the summary line the algorithms print |

For anything else geometric, use `compas.geometry` (`centroid_points`, `volume_polyhedron`, ...),
not a new helper.
