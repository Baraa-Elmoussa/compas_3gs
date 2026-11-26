# A lattice transmission tower with two cross-arms, under conductor and wind loads

The tower tapers from a 6 m square base to a 2 m square top over 24 m, in six panels: four legs,
horizontal struts at every level and X-bracing on every face. Two cross-arms carry the conductors.
The members are steel angles: the legs larger than the bracing and the horizontals.

Every face is X-braced, so the tower is highly statically indeterminate: statics leaves one free
parameter per self-stress state. `solve` gives the elastic state for the member sections. The
footing reactions do not depend on that choice: by statics, the two footings of each side share
W / 4 +- M / (2 b) of the weight W and of the overturning moment M about the base, b the base width.

Select a member, in the scene tree or in the view, to see its force in the info panel.

![A lattice transmission tower with two cross-arms, under conductor and wind loads](../assets/images/examples/model_06_transmission_tower.png)

```python title="examples/model/06_transmission_tower.py"
--8<-- "examples/model/06_transmission_tower.py"
```
