# $Id: shell_render_synth.py,v 1.4 2004/01/15 11:00:53 cpbotha Exp $
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


vr = vtk.vtkImageReader()
vr.SetHeaderSize(0)
vr.SetFileDimensionality(3)
vr.SetFileName("cube.raw")
vr.SetDataScalarType(3) # VTK_UNSIGNED_CHAR
vr.SetDataExtent((0,63,0,63,0,63))
#vr.SetDataSpacing((0.5, 3, 0.5))
vr.SetDataSpacing((0.5, 0.5, 0.5))

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
splatmapper.SetInput(vr.GetOutput())
splatmapper.SetRenderMode(1)

vprop = vtk.vtkVolumeProperty()
vprop.SetScalarOpacity(otf)
vprop.SetColor(ctf);
vprop.ShadeOn()
vprop.SetAmbient(0.1)
vprop.SetDiffuse(0.7)
vprop.SetSpecular(0.2)
vprop.SetSpecularPower(10)


volume = vtk.vtkVolume()
volume.SetProperty(vprop)
volume.SetMapper(splatmapper)

ren = vtk.vtkRenderer()
ren.AddVolume(volume)
#ren.GetActiveCamera().ParallelProjectionOn()

renwin = vtk.vtkRenderWindow()
renwin.AddRenderer(ren)

rwi = vtk.vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)

rwi.AddObserver('CharEvent', ce_cb)

#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()

