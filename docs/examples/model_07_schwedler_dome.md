# A Schwedler dome: rings, meridians and one diagonal per panel, on a ring of pinned supports

Schwedler domes, first built in Berlin in the 1860s, are spherical caps of straight bars: here a
30 m span and 10 m rise, 12 meridians, four rings and a crown node. The nodes carry the self-weight
of the cladding by their share of the dome surface, and in a second load case snow on one half.

Under the self-weight the dome works as a membrane: the meridians are in compression all the way
down, while the hoop forces in the rings change sign from compression near the crown to tension near
the base. Membrane theory puts that change at 51.8 degrees from the crown for a spherical dome under
its own weight. The dome is statically indeterminate; `solve` gives the elastic state for equal
member sections. Snow on one half puts the rings partly in tension and partly in compression.

Select a member, in the scene tree or in the view, to see its force in the info panel.

![A Schwedler dome: rings, meridians and one diagonal per panel, on a ring of pinned supports](../assets/images/examples/model_07_schwedler_dome.png)

```python title="examples/model/07_schwedler_dome.py"
--8<-- "examples/model/07_schwedler_dome.py"
```
