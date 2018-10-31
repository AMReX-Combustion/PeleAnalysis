.. highlight:: bash


stream - Streamlines of plotfile vector
***************************************

Given a plotfile containing the components of a vector field (velocity, gradient, etc), and an MEF file
containing a collection of "seed" points, create "streamlines" eminating from the seed points that are
locally parallel to the gradient field.  The resulting strealines will be of a fixed length going
both up and down the gradient from the seed point.  Results will written in a custom plotfile-like data folder
based on AMReX binary MultiFab I/O routines for portability, however the files ARE NOT PLOTFILES, and cannot
be processed or visualized with standard plotfile tools. The resulting stream lines retain the connectivity
inferred by the input seed points, which are expressed as a triangulated surface.  Thus, the streamlines
will bound a triangular-prism shaped volume that discretizes the layer bounded by a triangulated surface
obtained by connecting together the "low" ends of the streamlines, and a second surface connecting the "high"
ends of the streamlines.  Typically, this structure is used to discretize a region in space bounded by
two values of a scalar field - such as the reaction zone in a premixed flame.


```
Usage:
```

Example:

