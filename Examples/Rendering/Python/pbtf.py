import vtk
from vtk.util.colors import *

cubeSource = vtk.vtkCubeSource()
cubeSource.SetBounds(-5, 5, -5, 5, -5, 5)

ee = vtk.vtkExtractEdges()
ee.SetInput(cubeSource.GetOutput())

cubeTube = vtk.vtkTubeFilter()
cubeTube.SetRadius(0.1)
cubeTube.SetNumberOfSides(20)
cubeTube.UseDefaultNormalOn()
cubeTube.SetDefaultNormal(.5, .5, .5)
cubeTube.SetInput(ee.GetOutput())

cubeMapper = vtk.vtkPolyDataMapper()
cubeMapper.SetInput(cubeTube.GetOutput())

cubeActor = vtk.vtkActor()
cubeActor.SetMapper(cubeMapper)
cubeActor.GetProperty().SetDiffuseColor(lamp_black)
cubeActor.GetProperty().SetSpecular(0.4)
cubeActor.GetProperty().SetSpecularPower(10)

cubeMapper2 = vtk.vtkPolyDataMapper()
cubeMapper2.SetInput(cubeSource.GetOutput())
cubeActor2 = vtk.vtkActor()
cubeActor2.SetMapper(cubeMapper2)
cubeActor2.GetProperty().SetDiffuseColor(banana)
#cubeActor2.GetProperty().SetOpacity(0.9)

sphereSource = vtk.vtkSphereSource()
sphereSource.SetPhiResolution(32)
sphereSource.SetThetaResolution(32)
sphereSource.SetRadius(0.5)

sphereMapper = vtk.vtkPolyDataMapper()
sphereMapper.SetInput(sphereSource.GetOutput())

sphereActor = vtk.vtkActor()
sphereActor.SetMapper(sphereMapper)
sphereActor.GetProperty().SetDiffuseColor(tomato)


planeSize = 10
planeSources = [vtk.vtkPlaneSource() for dummy in range(3)]

for planeSource in planeSources:
    planeSource.SetOrigin(planeSize,planeSize,0)
    

# orthogonal to z axis
planeSources[0].SetPoint1(-planeSize, planeSize, 0)
planeSources[0].SetPoint2(planeSize, -planeSize, 0)

# orthogonal to y axis
planeSources[1].SetPoint1(-planeSize,planeSize,0)
planeSources[1].SetPoint2(planeSize,planeSize,planeSize*2)

# orthogonal to x axis
planeSources[2].SetPoint1(planeSize, planeSize, planeSize*2)
planeSources[2].SetPoint2(planeSize, -planeSize, 0)



planeMappers = [vtk.vtkPolyDataMapper() for dummy in planeSources]
for planeMapper,planeSource in zip(planeMappers, planeSources):
    planeMapper.SetInput(planeSource.GetOutput())

    
aRenderer = vtk.vtkRenderer()
aRenderer.SetBackground(slate_grey)
aRenderer.AddActor(cubeActor)
aRenderer.AddActor(cubeActor2)
aRenderer.AddActor(sphereActor)

planeActors = [vtk.vtkActor() for dummy in planeMappers]
for planeActor,planeMapper in zip(planeActors,planeMappers):
    planeActor.SetMapper(planeMapper)
    #planeActor.GetProperty().SetDiffuseColor(banana)
    planeActor.GetProperty().SetOpacity(.2)
    planeActor.GetProperty().SetAmbientColor(0,1,0)
    planeActor.GetProperty().SetAmbient(1)
    aRenderer.AddActor(planeActor)

planes = [vtk.vtkPlane() for dummy in planeSources]
cutters = [vtk.vtkCutter() for dummy in planes]
strippers = [vtk.vtkStripper() for dummy in cutters]
cutTubes = [vtk.vtkTubeFilter() for dummy in cutters]
cutMappers = [vtk.vtkPolyDataMapper() for dummy in cutters]
cutActors = [vtk.vtkActor() for dummy in cutters]

for plane,cutter,planeSource,stripper,cutTube,cutMapper,cutActor in zip(
    planes,cutters,planeSources,strippers,cutTubes,cutMappers,cutActors):
    
    cutter.SetCutFunction(plane)
    cutter.SetInput(cubeSource.GetOutput())
    stripper.SetInput(cutter.GetOutput())

    cutTube.SetRadius(0.1)
    cutTube.SetNumberOfSides(20)
    cutTube.UseDefaultNormalOn()
    cutTube.SetDefaultNormal(.5,.5,.5)
    cutTube.SetInput(stripper.GetOutput())    
    
    cutMapper.SetInput(cutTube.GetOutput())
    cutActor.SetMapper(cutMapper)

    cutActor.GetProperty().SetDiffuseColor(lamp_black)
    cutActor.GetProperty().SetSpecular(0.4)
    cutActor.GetProperty().SetSpecularPower(10)

    
    aRenderer.AddActor(cutActor)
    



renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(aRenderer)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
iren.Initialize()

pw = vtk.vtkPointWidget()
pw.PlaceWidget(-planeSize*2, planeSize*2, -planeSize*2, planeSize*2, -planeSize*2, planeSize*2)
pw.SetPosition(planeSize,planeSize,0)
pw.SetInteractor(iren)
pw.AllOff()
pw.On()

def observerPW(widget, eventName):
    sphereActor.SetPosition(widget.GetPosition())
    for planeSource,cutter in zip(planeSources,cutters):
        planeSource.SetCenter(widget.GetPosition())
        cutter.GetCutFunction().SetNormal(planeSource.GetNormal())
        cutter.GetCutFunction().SetOrigin(planeSource.GetOrigin())

# init
observerPW(pw, '')
# connect
pw.AddObserver('InteractionEvent', observerPW)



renWin.Render()
iren.Start()

