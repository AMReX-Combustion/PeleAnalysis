import vtk
from vtk import *
import string
import sys

f = open(sys.argv[1],'r')
outfile = sys.argv[2]

l = f.readline().strip()
tokens = l.split(' ')
if tokens[0] != 'VARIABLES':
    print('File not in correct format')
    sys.exit(0)
varnames = tokens[2:]
Nscalars = len(varnames) - 3

l = f.readline().strip()
params = {}
tokens = l.split(' ')
for token in tokens:
    pair = token.split('=')
    if len(pair)==2:
        params[pair[0]] = pair[1]

if (params['F'] != 'FEPOINT') | (params['ET'] != 'TRIANGLE'):
    print('File must be in Tecplot "FEPOINT TRIANGLE" formal')
    sys.exit(0)
    
Nnodes = int(params['N'])
Nelts = int(params['E'])

Points = vtk.vtkPoints()
PointScalars = []
for i in range(Nscalars):
    PointScalars.append(vtk.vtkDoubleArray())
    PointScalars[-1].SetName(varnames[3+i]);

print("Reading %d nodes..." % Nnodes)
for i in range(Nnodes):
    line = f.readline().strip(' \n')
    d = line.split(' ')
    d = [ x for x in d if x is not '' ]
    Nscalars_this = len(d) - 3
    if Nscalars != Nscalars_this:
        print('Number of scalars inconsistent across points')
        sys.exit(0)
    id = Points.InsertNextPoint(float(d[0]),float(d[1]),float(d[2]))
    for i in range(Nscalars):
        PointScalars[i].InsertNextTuple([float(d[3+i])])
    

print("Done (found position and %d scalar(s))" % Nscalars) 

Triangles = vtk.vtkCellArray()
Triangle = vtk.vtkTriangle()

print("Reading %d elements..." % Nelts)
for i in range(Nelts):
    line = f.readline().strip()
    d = line.split(' ')
    Triangle.GetPointIds().SetId(0,int(d[0])-1)
    Triangle.GetPointIds().SetId(1,int(d[1])-1)
    Triangle.GetPointIds().SetId(2,int(d[2])-1)
    Triangles.InsertNextCell(Triangle)

print("Done") 
f.close()


polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
polydata.SetPolys(Triangles)
if (Nscalars>0):
    for scalar in PointScalars:
        polydata.GetPointData().AddArray(scalar)
    polydata.GetPointData().SetActiveScalars(varnames[-1])

polydata.Modified()
if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()
 
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName(outfile);
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.SetDataModeToBinary()
#writer.SetDataModeToAscii()
print("Writing data to %s ..." % outfile)
writer.Write()
