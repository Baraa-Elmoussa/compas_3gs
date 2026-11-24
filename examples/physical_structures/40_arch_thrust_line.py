"""Thrust line of a semicircular masonry arch, and what 3GS can and cannot do with it.

1. The arch: a semicircle of voussoirs (unit depth); each voussoir weighs its area times DENSITY.
2. 3GS first: at every point of a thrust line, three forces meet (the thrust from the left, the
   thrust to the right, the weight). Three forces in equilibrium always lie in one plane, and a
   polyhedral cell needs forces that span 3D. So 3GS has no force polyhedron for an arch node:
   `egi_from_vectors` refuses fewer than four vectors, and four or more forces in one plane do not
   close a polyhedron either. The thrust line of a planar arch is 2D graphic statics.
3. The thrust line with 2D graphic statics: the funicular polygon of the voussoir weights for a
   horizontal thrust H, through points of the two springing joints at x = -d and +d. The arch is
   safe (Heyman) for every thrust line that stays inside the arch at every joint. Scanning H and d
   gives the range between the minimum and the maximum thrust.
4. Viewer: the arch, the thrust lines of the minimum and maximum thrust, and the force polygon
   (load line and rays from the pole) of the minimum-thrust state.

For thrust lines of 3D masonry (vaults, domes) the BRG tool is compas_tna (thrust network analysis);
3GS is for spatial structures where nodes have four or more forces that do not lie in one plane.

Run with the compas-3gs env:  python examples/physical_structures/40_arch_thrust_line.py

"""

import math

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Polygon
from compas.geometry import Polyline
from compas.geometry import intersection_segment_segment
from compas_3gs.algorithms import egi_from_vectors

RADIUS = 5.0  # centreline radius
THICKNESS = 1.0
VOUSSOIRS = 15
DENSITY = 20.0  # kN/m3, so weights are kN per metre of depth
H_STEP = 0.5  # kN/m, resolution of the scan over the horizontal thrust
D_STEPS = 20  # positions in each springing joint

# ------------------------------------------------------------------------------
# 1. voussoirs, joints and weights
# ------------------------------------------------------------------------------

ri, re = RADIUS - THICKNESS / 2, RADIUS + THICKNESS / 2
angles = [math.pi - k * math.pi / VOUSSOIRS for k in range(VOUSSOIRS + 1)]  # left springing to right
joints = [([ri * math.cos(a), ri * math.sin(a), 0], [re * math.cos(a), re * math.sin(a), 0]) for a in angles]

weights = []
xs = []
for k in range(VOUSSOIRS):
    half = abs(angles[k] - angles[k + 1]) / 2
    area = half * (re**2 - ri**2)
    rc = 2 / 3 * (re**3 - ri**3) / (re**2 - ri**2) * math.sin(half) / half  # centroid radius of an annular sector
    weights.append(area * DENSITY)
    xs.append(rc * math.cos((angles[k] + angles[k + 1]) / 2))
W = sum(weights)
print("arch: radius {}, thickness {} (t/R = {:.3f}), {} voussoirs, total weight {:.1f} kN/m".format(RADIUS, THICKNESS, THICKNESS / RADIUS, VOUSSOIRS, W))

# ------------------------------------------------------------------------------
# 2. 3GS at a node of the thrust line
# ------------------------------------------------------------------------------

H_test = 50.0
node_forces = {0: [H_test, W / 2, 0.0], 1: [0.0, -weights[0], 0.0], 2: [-H_test, -(W / 2 - weights[0]), 0.0]}
try:
    egi_from_vectors(node_forces, [0, 0, 0])
except Exception as error:
    print("3GS with the 3 forces at a thrust-line node: {}".format(error))
node_forces[3] = [0.0, 0.0, 1e-4]  # a fourth, tiny force out of the plane does not help
try:
    egi_from_vectors(node_forces, [0, 0, 0])
except Exception as error:
    print("3GS with a 4th force added: {}".format(error))

# ------------------------------------------------------------------------------
# 3. thrust lines (funicular polygons) for a range of H and springing positions d
# ------------------------------------------------------------------------------

admissible = {}  # (H, d): thrust-line points, for the thrust lines inside the arch
for i in range(D_STEPS + 1):
    d = ri + i * (re - ri) / D_STEPS
    for j in range(int(10 / H_STEP), int(200 / H_STEP)):
        H = j * H_STEP
        # funicular of the weights from (-d, 0) to (+d, 0): left reaction from moments about +d
        shear = sum(w * (d - x) for w, x in zip(weights, xs)) / (2 * d)
        points = [[-d, 0.0, 0.0]]
        for w, x in zip(weights, xs):
            points.append([x, points[-1][1] + shear / H * (x - points[-1][0]), 0.0])
            shear -= w
        points.append([d, points[-1][1] + shear / H * (d - points[-1][0]), 0.0])
        # inside the arch if it crosses every joint between its intrados and extrados points
        segments = list(zip(points[:-1], points[1:]))
        if all(any(intersection_segment_segment(joint, segment, tol=1e-9)[0] for segment in segments) for joint in joints[1:-1]):
            admissible[(H, d)] = points

if admissible:
    centred = [H for H, d in admissible if abs(d - RADIUS) < 1e-9]
    (H_min, d_min), (H_max, d_max) = min(admissible), max(admissible)
    if centred:
        print("thrust through the middle of both springing joints: H from {:.1f} to {:.1f} kN/m".format(min(centred), max(centred)))
    print("thrust anywhere in the springing joints: H from {:.1f} kN/m (x = ±{:.2f}) to {:.1f} kN/m (x = ±{:.2f})".format(H_min, d_min, H_max, d_max))
    print("  minimum thrust H/W = {:.3f}, maximum thrust H/W = {:.3f}".format(H_min / W, H_max / W))
else:
    print("no thrust line fits inside the arch: it cannot stand (t/R too small)")

# ------------------------------------------------------------------------------
# 4. viewer: arch, thrust lines, force polygon
# ------------------------------------------------------------------------------

viewer = Viewer()
for k in range(VOUSSOIRS):
    (a, b), (c, d) = joints[k], joints[k + 1]
    viewer.scene.add(Polygon([a, b, d, c]), facecolor=Color.grey().lightened(60), linecolor=Color.grey())

if admissible:
    for key, color, name in (((H_min, d_min), Color.red(), "H min"), ((H_max, d_max), Color.blue(), "H max")):
        points = admissible[key]
        viewer.scene.add(Polyline(points), linecolor=color, linewidth=3)
        viewer.scene.add(Tag("{} = {:.1f} kN/m".format(name, key[0]), [0, max(p[1] for p in points) + 0.6, 0], color=color, height=16))

    # force polygon of the minimum-thrust state, to the right of the arch (1 m drawn per 10 kN)
    top = [RADIUS + 3 + H_min / 10, W / 20, 0]
    pole = [top[0] - H_min / 10, top[1] - sum(w * (d_min - x) for w, x in zip(weights, xs)) / (2 * d_min) / 10, 0]
    viewer.scene.add(Point(*pole), pointcolor=Color.red(), pointsize=10)
    viewer.scene.add(Line(pole, top), linecolor=Color.red())
    y = top[1]
    for w in weights:
        viewer.scene.add(Line([top[0], y, 0], [top[0], y - w / 10, 0]), linecolor=Color.black(), linewidth=3)
        y -= w / 10
        viewer.scene.add(Line(pole, [top[0], y, 0]), linecolor=Color.red())
    viewer.scene.add(Tag("force polygon (1 m = 10 kN)", [top[0] - 1, top[1] + 0.8, 0], color=Color.black(), height=14))

viewer.show()
