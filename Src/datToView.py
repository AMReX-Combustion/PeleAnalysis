import vtk
from vtk import *
import string
import sys

def ReadZoneHeader(f):
    l = f.readline().strip()
    params = {}
    tokens = l.split(' ')
    for token in tokens:
        pair = token.split('=')
        if len(pair)==2:
            params[pair[0]] = pair[1].strip('\'')

    NnodesPerElt = -1
    try:
        if params['F'] == 'FEPOINT':
            try:
                if params['ET'] == 'TRIANGLE':
                    #print('Found ET=TRIANGLE F=FEPOINT')
                    NnodesPerElt = 3
                elif params['ET'] == 'SEGMENT':
                    #print('Found ET=SEGMENT F=FEPOINT')
                    NnodesPerElt = 2
                else:
                    print('F=FEPOINT requires ET=TRIANGLE format specifier')
                    sys.exit(0)
            except KeyError:
                print('F=FEPOINT requires ET=TRIANGLE format specifier')
                sys.exit(0)
    except KeyError:
        try:
            if (params['DATAPACKING'] == 'POINT'):
                try:
                    if (params['ZONETYPE'] == 'FETRIANGLE'):
                        #print('Found ZONETYPE=FETRIANGLE DATAPACKING=POINT')
                        NnodesPerElt = 3
                    elif (params['ZONETYPE'] == 'FELINESEG'):
                        #print('Found ZONETYPE=FELINESEG DATAPACKING=POINT')
                        NnodesPerElt = 2
                    else:
                        print('DATAPACKING=POINT requires ZONETYPE=FETRIANGLE format specifier')
                except KeyError:
                    print('DATAPACKING=POINT requires ZONETYPE=FETRIANGLEformat specifier')
            else:            
                print('Neither F=FEPOINT nor DATAPACKING=POINT format specified in input')
        except KeyError:
            return 0, 0, 0
            #print('Neither F or DATAPACKING format specified in input')

    Nnodes = int(params['N'])
    Nelts = int(params['E'])
    return Nnodes, Nelts,NnodesPerElt

def ReadNodes(Nnodes,Nscalars,varnames,f):
    Points = vtk.vtkPoints()
    PointScalars = []
    for i in range(Nscalars):
        PointScalars.append(vtk.vtkDoubleArray())
        PointScalars[-1].SetName(varnames[3+i]);

    for i in range(Nnodes):
        line = f.readline().strip(' \n')
        d = line.split(' ')
        d = [ x for x in d if x is not '' ]
        Nscalars_this = len(d) - 3
        if Nscalars != Nscalars_this:
            print('Number of scalars inconsistent across points '+str(Nscalars)+' vs '+str(Nscalars_this))
            print(d)
            sys.exit(0)
        id = Points.InsertNextPoint(float(d[0]),float(d[1]),float(d[2]))
        for i in range(Nscalars):
            PointScalars[i].InsertNextTuple([float(d[3+i])])

    return Points, PointScalars
    
def ReadElts(Nelts,NnodesPerElt,f):
    polys = vtk.vtkCellArray()

    if NnodesPerElt==2:
        #print('Reading segments')
        for i in range(Nelts):
            polys.InsertNextCell(2)
            polys.InsertCellPoint(i)
            polys.InsertCellPoint(i+1)
            line = f.readline().strip() # Read and toss, do not need for PolyLine structs

    elif NnodesPerElt==3:
        #print('Reading triangles')
        Triangle = vtk.vtkTriangle()
        for i in range(Nelts):
            line = f.readline().strip()
            d = line.split(' ')
            Triangle.GetPointIds().SetId(0,int(d[0])-1)
            Triangle.GetPointIds().SetId(1,int(d[1])-1)
            Triangle.GetPointIds().SetId(2,int(d[2])-1)
            polys.InsertNextCell(Triangle)

    return polys
 

f = open(sys.argv[1],'r')

l = f.readline().strip()
tokens = l.split(' ')
if tokens[0] != 'VARIABLES':
    print('File not in correct format')
    sys.exit(0)
varnames = tokens[2:]
Nscalars = len(varnames) - 3

allpolys = vtk.vtkAppendPolyData()
EOF = False
while(not EOF):
    try:
        Nnodes, Nelts, NnodesPerElt = ReadZoneHeader(f)
        if Nnodes==0 | Nelts==0:
            EOF = True
        else:
            #print('has Nnodes, Nelts: '+str(Nnodes)+' '+str(Nelts))
            Points, PointScalars = ReadNodes(Nnodes, Nscalars,varnames,f)
            polys = ReadElts(Nelts,NnodesPerElt,f)
            
            polydata = vtk.vtkPolyData()
            polydata.SetPoints(Points)
            if NnodesPerElt==2:
                polydata.SetLines(polys)
            else:
                polydata.SetPolys(polys)
            if (Nscalars>0):
                for scalar in PointScalars:
                    polydata.GetPointData().AddArray(scalar)
                polydata.GetPointData().SetActiveScalars(varnames[-1])
                    
            polydata.Modified()
            if vtk.VTK_MAJOR_VERSION <= 5:
                polydata.Update()

            allpolys.AddInputData(polydata)
            allpolys.Update()

    except:
       EOF=True

mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(allpolys)
