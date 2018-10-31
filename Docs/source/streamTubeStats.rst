.. highlight:: bash


streamTubeStats
***************

Given a set of streamlines that represents the bounds of triangular-wedge shaped volumes, gather statistics
of the volumes.  Output the results as a triangulated surface MEF file, where for each triangle, the values at all three
nodes are equal and represent the quanity computed for the entire wedge volume.  Thus, the MEF format is overloaded
here so that each node is multiply-defined, based on the number of triangles is it part of. The resulting MEF
file will also contain the area of the triangle from the original tesselation that created the seed points
for the streamline, with the intention that the results here can be used to construct statistics weighted by
this area.  Typically, this processing is used to gather statistics associated with a isosurface (such as one that
represents a flame).


```
Usage:
```

Example:

