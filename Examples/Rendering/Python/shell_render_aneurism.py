# $Id: shell_render_aneurism.py,v 1.4 2004/01/15 11:00:53 cpbotha Exp $
# example to test shell renderer (*shudder*)

from vtkpython import *
from vtkdevide import *
import time

def bench(camera, rwi):
    initial_time = time.clock()
    for i in range(36):
        camera.Azimuth(10)
        rwi.Render()
    
    end_time = time.clock()

    print "FPS == %f" % (36 / (end_time - initial_time))

def ce_cb(obj, evt_name):
    if obj.GetKeyCode() == 'm':
        crm = splatmapper.GetRenderMode()
	crm = crm + 1
	if crm > 2:
	    crm = 0
        splatmapper.SetRenderMode(crm)
        print "rendermode switched to %d" % (crm)

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
        bench(ren.GetActiveCamera(), rwi)
        
    rwi.Render()
        
        
reader = vtkImageReader()
reader.SetHeaderSize(0)
reader.SetFileDimensionality(3)
reader.SetFileName("aneurism.raw")
reader.SetDataScalarType(3)
reader.SetDataExtent(0,255,0,255,0,255)
reader.SetDataSpacing(1,1,1)

otf = vtkPiecewiseFunction()
otf.AddPoint(0.0, 0.0)
otf.AddPoint(100, 0.0)
otf.AddPoint(100.1, 0.9)
otf.AddPoint(255.0, 0.9)

ctf = vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(100, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(100.1, 0.82, 0.43, 0.35)
ctf.AddRGBPoint(255.0, 0.82, 0.43, 0.35)

#se = vtkShellExtractor()
#se.SetInput(hdfr.GetOutput())
#se.SetOpacityTF(otf)
#se.SetOmegaL(0.8)
#se.SetOmegaH(0.99)

#se.Update()

splatmapper = vtkOpenGLVolumeShellSplatMapper()
splatmapper.SetOmegaL(0.1)
splatmapper.SetOmegaH(0.2)
splatmapper.SetInput(reader.GetOutput())
splatmapper.SetRenderMode(0)

vprop = vtkVolumeProperty()
vprop.SetScalarOpacity(otf)
vprop.SetColor(ctf);
vprop.ShadeOn()
vprop.SetAmbient(0.1)
vprop.SetDiffuse(0.7)
vprop.SetSpecular(0.2)
vprop.SetSpecularPower(10)

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
cubeAxesActor2d.VisibilityOn()
reader.Update()
cubeAxesActor2d.SetBounds(reader.GetOutput().GetBounds())
cubeAxesActor2d.SetCamera(ren.GetActiveCamera())


renwin = vtkRenderWindow()
renwin.AddRenderer(ren)

rwi = vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)

rwi.AddObserver('CharEvent', ce_cb)

#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()

