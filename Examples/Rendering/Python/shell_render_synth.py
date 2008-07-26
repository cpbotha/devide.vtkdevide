# $Id$
# example to test shell renderer (*shudder*)

import vtk
import vtkdevide

def ce_cb(obj, evt_name):
    if obj.GetKeyCode() == 'm':
        crm = splatmapper.GetRenderMode()
        splatmapper.SetRenderMode(not crm)
        print "rendermode switched to %d" % (not crm)
        
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
        
    rwi.Render()


reader = vtk.vtkXMLImageDataReader()
reader.SetFileName("cube.vti")

otf = vtk.vtkPiecewiseFunction()
otf.AddPoint(0.0, 0.0)
otf.AddPoint(253.9, 0.0)
otf.AddPoint(254, 0.95)
otf.AddPoint(255, 0.95)

ctf = vtk.vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(253.9, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(254, 1.0, 1.0, 1.0)
ctf.AddRGBPoint(254.9, 1.0, 1.0, 1.0)
ctf.AddRGBPoint(255, 1.0, 0, 0)

#se = vtkShellExtractor()
#se.SetInput(hdfr.GetOutput())
#se.SetOpacityTF(otf)
#se.SetOmegaL(0.8)
#se.SetOmegaH(0.99)

#se.Update()

splatmapper = vtkdevide.vtkOpenGLVolumeShellSplatMapper()
splatmapper.SetOmegaL(0.9)
splatmapper.SetOmegaH(0.9)
splatmapper.SetInput(reader.GetOutput())
splatmapper.SetRenderMode(1)
splatmapper.SetPerspectiveOrderingMode(1)

vprop = vtk.vtkVolumeProperty()
vprop.SetScalarOpacity(otf)
vprop.SetColor(ctf);
vprop.ShadeOn()
vprop.SetAmbient(0.4)
vprop.SetDiffuse(0.7)
vprop.SetSpecular(0.2)
vprop.SetSpecularPower(70)


volume = vtk.vtkVolume()
volume.SetProperty(vprop)
volume.SetMapper(splatmapper)

ren = vtk.vtkRenderer()
ren.AddVolume(volume)
#ren.GetActiveCamera().ParallelProjectionOn()

cubeAxesActor2d = vtk.vtkCubeAxesActor2D()
cubeAxesActor2d.SetFlyModeToOuterEdges()
ren.AddActor(cubeAxesActor2d)
cubeAxesActor2d.VisibilityOn() # FIXME: activate axes here
reader.Update()
cubeAxesActor2d.SetBounds(reader.GetOutput().GetBounds())
cubeAxesActor2d.SetCamera(ren.GetActiveCamera())

renwin = vtk.vtkRenderWindow()
renwin.AddRenderer(ren)

rwi = vtk.vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)

rwi.AddObserver('CharEvent', ce_cb)

#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()

