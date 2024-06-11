"""
Converts Isosurface data in Marc's Element Format
to the VTK format.
The resulting file is a bit larger than MEF but the conversion
takes only a few seconds even for meshes with > 30M points.
"""
import sys
import numpy as np
import meshio

def shape_from_header(h):
    """
    Infer the shape the FAB data and the number of fields
    from the header in a plotfile binary file
    h: string of the header line
    """
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    total_shape = [shape[0], shape[1], shape[2], nfields]
    return total_shape


mef_file = sys.argv[1]

with open(mef_file, 'rb') as mef:
    # Read Marc's DNA
    time = mef.readline().decode('ascii')
    fields = mef.readline().decode('ascii').split()
    nfaces, verts_per_face = mef.readline().decode('ascii').split()
    nfaces = int(nfaces)
    verts_per_face = int(verts_per_face)
    # Read and parse the FAB header
    header = mef.readline().decode('ascii')
    bin_data_shape = shape_from_header(header)
    #  Handle the 1D format of the data
    bin_data_shape = (bin_data_shape[0], bin_data_shape[3])
    # This can take time
    print('Reading data...')
    # Read the n_verts x n_fields floats
    data = np.fromfile(mef, 'float64', np.prod(bin_data_shape))
    # Not Fortran ordered like AMReX
    data = data.reshape(bin_data_shape, order='C')
    # Read the nfaces x verts_per_face ints from bytes
    int_data_shape = (nfaces, verts_per_face)
    faces = np.fromfile(mef, 'int32', np.prod(int_data_shape))
    # Reshape into array of faces (indexed at 0)
    faces = faces.reshape(int_data_shape, order='C') - 1

# This is quite fast
print('Creating mesh...')
# Infer the number of dimensions from the field names
dim_fields = ['X', 'Y', 'Z']
ndims = 0
for f in fields:
    if f in dim_fields:
        ndims += 1
# Array of vertices
vertices = data[:, :ndims]
# Dictionnary of point data to have field
# names in the .vtk file
point_data = data[:, ndims:]
point_data = {f:point_data[:, i] for i, f in enumerate(fields[ndims:])}
# Create the mesh object
mesh = meshio.Mesh(points=vertices,
                   # All faces are triangles
                   cells={'triangle':faces},
                   point_data=point_data)

# This is also quite fast
print('Saving mesh object...')
# Just change the extension to '.vtk'
mesh.write('.'.join(mef_file.split('.')[:-1] + ['vtk']))




