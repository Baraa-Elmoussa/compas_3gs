# Interpreting diagrams

A force **f**~i,j~ acting on node vi can be either a compression force or a tension force; it is a compression force if it is pushing onto the node, or a tension force if it is pulling away from the node.
In 3D graphic statics, the direction of **f**~i,j~ is the same as the direction of the normal **n̂** ^\*^~i,j~ of the corresponding face of the polyhedral cell.
Therefore, a single force polyhedron can represent the equilibrium of various combinations of compression and tension forces ([^1]).
For example in b1 and b2 in the figure below, the directions of the external forces are in the same orientations and directions as the corresponding face normals in *c*^\*^~i~ (a in figure below).
However, the external forces are interpreted differently; **f**~i,4~ of $\Lambda$ case 1 is in tension, while **f**~i,4~ of $\Lambda$ case 2 is in compression.

The interpretation the *j*-th external force or the internal force of the *j*-th
structural member at node **v**~i~ as either compression or tension, can be made by comparing the face normal **n̂** ^\*^~i,j~ and the orientation of the corresponding member ([^2]).
Locally defined for each node **v**~i~, the orientation of the *j*-th external force or member can be represented by a unitised edge vector **\`ê\`**~i,j~ with the tail of the vector at **v**~i~ (c in figure below).
If **n̂** ^\*^~i,j~ and **ê**~i,j~ are pointing in the same direction, the force in the corresponding member is interpreted as positive, so in tension; if the vectors are in opposite directions, the force in the corresponding member is interpreted as negative, so in compression.
The compression and tension forces are shown in blue and red, respectively.

<br>


![](../../assets/images/cell_force_interpretation.jpg){width="100.0%"}
/// caption
Multiple configurations of compression and tension elements are possible for a polyhedral cell: a) $\Lambda^{\perp}$ or polyhedral cell *c*^\*^~i~ of node **v**~i~; b) orientation of the members or external forces as unitised vectors **ê**~i,j~; and c) the direction and type of forces in the members.
///

<br>


------------------------------------------------------------------------

## References

[^1]: Akbarzadeh, M. (2016). *3D Graphical Statics using 3D Reciprocal Diagrams*. Ph.D. thesis, ETH Zürich, Zürich, Switzerland.

[^2]: Lee, J., T. Van Mele, and P. Block (2016). Form-finding explorations through geometric transformations and modifications of force polyhedrons. In *Proceedings of the 2016 International Association for Shell and Spatial Structures (IASS) Symposium*, Tokyo, Japan.
