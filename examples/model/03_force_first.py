"""Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces.

A cube split into six pyramids is the force diagram of an octahedral frame under six loads.
`GSModel.from_force_diagram` makes the dual network the form diagram, reciprocates it so that every
member is perpendicular to its face, takes the member forces from the face areas and the loads
from the boundary faces.

Run with the compas-3gs env:  python examples/model/03_force_first.py

"""

from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.model import GSModel
from compas_3gs.operations import volmesh_cell_subdivide_barycentric
from compas_3gs.viewer import GSViewer

force = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
volmesh_cell_subdivide_barycentric(force, next(iter(force.cells())))

model = GSModel.from_force_diagram(force)

print(model.summary())
for node in model.form.nodes():
    print("node {} load {}".format(node, [round(component, 3) for component in model.form.node_load(node)]))

viewer = GSViewer()
viewer.add_model(model)
viewer.show()
