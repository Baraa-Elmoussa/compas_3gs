# A tensegrity prism: three struts held apart by nine cables, in equilibrium without any load

Kenneth Snelson's simplest tensegrity, the module of his Needle Tower: two equilateral triangles of
cables, the top one turned against the bottom one, three vertical cables and three struts.

With no load, equilibrium can only be a self-stress state: a prestress that the struts and cables
hold in balance among themselves. Statics finds it: the prism has one self-stress state, and only
when the top triangle is turned by 150 degrees; at any other angle it has none and cannot be
prestressed. The level of the prestress is free, so it is prescribed through one strut force. Then
all struts are in compression and all cables in tension, as a tensegrity needs.

Every node has four forces, one strut and three cables, so its force polyhedron is a tetrahedron.

The prism is also a mechanism to first order: it carries an external load through the stiffness its
prestress gives it, which is beyond the linear statics of this example.

![A tensegrity prism: three struts held apart by nine cables, in equilibrium without any load](../assets/images/examples/model_08_tensegrity_prism.png)

```python title="examples/model/08_tensegrity_prism.py"
--8<-- "examples/model/08_tensegrity_prism.py"
```
