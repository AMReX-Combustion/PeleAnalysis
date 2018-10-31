.. highlight:: bash


surfMEF Tools
*************

The MEF (Marc's element format) file format is hacked up data structure-on-disk, primarily intended to store
a triangulated surface.  The data at the nodes (including the position, but also potentially other field
quantities) is written using AMReX's FabIO functions so that the floating point data is portable.  The
triangles are specified as a set of integer triples, where each integer is the (zero-based) id of the
node.  The first four lines in an MEF file are ASCII, adn include a label, the list of variables,
the integer number of elements, and then the FAB header info.  The binary FAB info is then concatenated, and
is then followed by the element triples, one element per one, written in ASCII.

The surfMEF conversion tools are used to transform to and from the MEF format, and to do simple arithmetic
operations on the data.

- combineMEF: Combine the components of two MEF files, assuming they have the same node positions and connectivity
- mergeMEF: Merge the triangles of two different MEF files
- multMEF: Multiply specific components of MEF files together
- scaleMEF: Scale specific components of the MEF by constants
- sliceMEF: Compute a contour on an MEF surface
- smoothMEF: Smooth an MEF surface
- surfDATtoMEF: Convert a Tecplot-formatted ASCII triangulated surface file into an MEF file
- surfMEFtoDAT: Convert an MEF-file into a Tecplot-formatted ASCII triangulated surface file
- trimMEFgen: Do an area-weighted binning of an MEF surface file, by assuming linear variatoin of the field on each triagle and slicing the triangles into bits at the bin boundaries


```
Usage:
```

Example:

