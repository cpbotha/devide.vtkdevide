# $Id: shell_render_stanford_head.py,v 1.3 2004/01/14 18:51:38 cpbotha Exp $
# example to test shell renderer (*shudder*)

from vtk import *
from vtkdscas import *
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
        bench(ren.GetActiveCamera(), rwi)
        
	
        
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
otf.AddPoint(1000.9, 0.0)
otf.AddPoint(1001, 0.2)
otf.AddPoint(1249.9, 0.2)
otf.AddPoint(1250.0, 1)
otf.AddPoint(2800.0, 1)

ctf = vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(1000.9, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(1001, 0.0, 1.0, 0.0)
ctf.AddRGBPoint(1249.9, 0.0, 1.0, 0.0)
ctf.AddRGBPoint(1250, 1.0, 0.937, 0.859)
ctf.AddRGBPoint(2800, 1.0, 0.937, 0.859)

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
vprop.SetSpecular(0.2)
vprop.SetSpecularPower(10)

volume = vtkVolume()
volume.SetProperty(vprop)
volume.SetMapper(splatmapper)

ren = vtkRenderer()
ren.SetBackground(0.5, 0.5, 0.5)
ren.AddVolume(volume)
#ren.GetActiveCamera().ParallelProjectionOn()

renwin = vtkRenderWindow()
renwin.AddRenderer(ren)

rwi = vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)

rwi.AddObserver('CharEvent', ce_cb)

#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()

