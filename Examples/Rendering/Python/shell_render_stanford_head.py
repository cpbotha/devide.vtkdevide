# $Id: shell_render_stanford_head.py,v 1.6 2004/06/22 14:26:48 cpbotha Exp $
# example to test shell renderer (*shudder*)

from vtk import *
import vtk
from vtkdevide import *
import time

def bench(camera, rwi):
    initial_time = time.clock()
    for i in range(36):
        camera.Azimuth(10)
        rwi.Render()
    
    end_time = time.clock()

    print "FPS == %f" % (36 / (end_time - initial_time))
    
def bench2(camera, rwi):
    initial_time = time.clock()

    numberOfRenders = 10 * (36 + 1)
    
    for i in range(10):
        for j in range(36):
            camera.Azimuth(10)
            rwi.Render()
        
        camera.Elevation(36 * i)
        rwi.Render()
            
            
    
    end_time = time.clock()

    print "FPS == %f" % (numberOfRenders / (end_time - initial_time))
    

def ce_cb(obj, evt_name):
    if obj.GetKeyCode() == 'm':
        crm = splatmapper.GetRenderMode()
        splatmapper.SetRenderMode(not crm)
        print "rendermode switched to %d" % (not crm)
        
    if obj.GetKeyCode() == 'i':
        com = splatmapper.GetPerspectiveOrderingMode()
	com = com + 1
	if com > 2:
	    com = 0
        splatmapper.SetPerspectiveOrderingMode(com)
        print "ordering mode switched to %d" % (com)


    elif obj.GetKeyCode() == '\'':
        cur = splatmapper.GetEllipsoidDiameter()
        splatmapper.SetEllipsoidDiameter(cur - 0.1)
        print "EllipsoidDiameter == %s" % str(cur - 0.1)
    elif obj.GetKeyCode() == ',':
        cur = splatmapper.GetEllipsoidDiameter()
        splatmapper.SetEllipsoidDiameter(cur + 0.1)
        print "EllipsoidDiameter == %s" % str(cur + 0.1)
        
    elif obj.GetKeyCode() == 'd':
        cur = splatmapper.GetGaussianRadialExtent()
        splatmapper.SetGaussianRadialExtent(cur - 0.1)
        print "GaussianRadialExtent == %s" % str(cur - 0.1)
    elif obj.GetKeyCode() == 'h':
        cur = splatmapper.GetGaussianRadialExtent()
        splatmapper.SetGaussianRadialExtent(cur + 0.1)
        print "GaussianRadialExtent == %s" % str(cur + 0.1)
        
    elif obj.GetKeyCode() == 't':
        cur = splatmapper.GetGaussianSigma()
        splatmapper.SetGaussianSigma(cur - 0.1)
        print "GaussianSigma == %s" % str(cur - 0.1)
    elif obj.GetKeyCode() == 'n':
        cur = splatmapper.GetGaussianSigma()
        splatmapper.SetGaussianSigma(cur + 0.1)
        print "GaussianSigma == %s" % str(cur + 0.1)
	
    elif obj.GetKeyCode() == 'b':
        bench2(ren.GetActiveCamera(), rwi)
        
	
        
    rwi.Render()
        
        
        

v16 = vtkVolume16Reader()
v16.SetDataDimensions(256,256)
v16.SetDataByteOrderToBigEndian()
v16.SetFilePrefix("CThead/CThead")
v16.SetImageRange(1, 99)
v16.SetDataSpacing(1, 1, 2)
v16.Update()

otf = vtkPiecewiseFunction()
otf.AddPoint(0.0, 0.0)
otf.AddPoint(899.9, 0.0)
otf.AddPoint(900, 0)
otf.AddPoint(1499.9, 1)
otf.AddPoint(1500.0, 1)
otf.AddPoint(65535.0, 1)

#skinCol = (0.83, 0.64, 0.58)
skinCol = (0.93, 0.87, 0.80)
boneCol = skinCol
#boneCol = (1.0, 0.937, 0.859)

ctf = vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(899.9, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(900, skinCol[0], skinCol[1], skinCol[2])
ctf.AddRGBPoint(1499.9, skinCol[0], skinCol[1], skinCol[2])
ctf.AddRGBPoint(1500, boneCol[0], boneCol[1], boneCol[2])
ctf.AddRGBPoint(2800, boneCol[0], boneCol[1], boneCol[2])

#se = vtkShellExtractor()
#se.SetInput(hdfr.GetOutput())
#se.SetOpacityTF(otf)
#se.SetOmegaL(0.8)
#se.SetOmegaH(0.99)

#se.Update()

splatmapper = vtkOpenGLVolumeShellSplatMapper()
splatmapper.SetOmegaL(0.9)
splatmapper.SetOmegaH(0.9)
splatmapper.SetInput(v16.GetOutput())
splatmapper.SetRenderMode(0)

#splatmapper = vtkVolumeTextureMapper2D()
#splatmapper.SetInput(v16.GetOutput())

vprop = vtkVolumeProperty()
vprop.SetScalarOpacity(otf)
vprop.SetColor(ctf);
vprop.ShadeOn()
vprop.SetAmbient(0.1)
vprop.SetDiffuse(0.7)
vprop.SetSpecular(0.4)
vprop.SetSpecularPower(60) # 10

volume = vtkVolume()
volume.SetProperty(vprop)
volume.SetMapper(splatmapper)

ren = vtkRenderer()
ren.SetBackground(0.5, 0.5, 0.5)
ren.AddVolume(volume)
#ren.GetActiveCamera().ParallelProjectionOn()

cubeAxesActor2d = vtk.vtkCubeAxesActor2D()
cubeAxesActor2d.SetFlyModeToOuterEdges()
ren.AddActor(cubeAxesActor2d)
cubeAxesActor2d.VisibilityOff() # turn on here!
v16.Update()
cubeAxesActor2d.SetBounds(v16.GetOutput().GetBounds())
cubeAxesActor2d.SetCamera(ren.GetActiveCamera())

renwin = vtkRenderWindow()
renwin.SetSize(512, 512)
renwin.AddRenderer(ren)

rwi = vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)

rwi.AddObserver('CharEvent', ce_cb)

#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()

