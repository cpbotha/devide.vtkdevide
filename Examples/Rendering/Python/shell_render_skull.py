# $Id: shell_render_skull.py,v 1.1 2003/01/08 14:07:29 cpbotha Exp $
# example to test shell renderer (*shudder*)

from vtkpython import *
from vtkcpbothapython import *

hdfr = vtkHDFVolumeReader()
hdfr.SetFileName("skull_256x256x256.hdf")

otf = vtkPiecewiseFunction()
otf.AddPoint(0.0, 0.0)
otf.AddPoint(39.9, 0.0)
otf.AddPoint(40.0, 1.0)
otf.AddPoint(150.0, 1.0)

ctf = vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(39.9, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(40, 1.0, 0.937, 0.859)
ctf.AddRGBPoint(150, 1.0, 0.937, 0.859)

#se = vtkShellExtractor()
#se.SetInput(hdfr.GetOutput())
#se.SetOpacityTF(otf)
#se.SetOmegaL(0.8)
#se.SetOmegaH(0.99)

#se.Update()

splatmapper = vtkOpenGLVolumeShellSplatMapper()
splatmapper.SetOmegaL(0.9)
splatmapper.SetOmegaH(0.9)
splatmapper.SetInput(hdfr.GetOutput())
splatmapper.SetRenderMode(0)

vprop = vtkVolumeProperty()
vprop.SetScalarOpacity(otf)
vprop.SetColor(ctf);
vprop.ShadeOn()
vprop.SetAmbient(0.4)
vprop.SetDiffuse(0.7)
#vprop.SetSpecular(0.2)
vprop.SetSpecularPower(70)

volume = vtkVolume()
volume.SetProperty(vprop)
volume.SetMapper(splatmapper)

ren = vtkRenderer()
ren.AddVolume(volume)
ren.GetActiveCamera().ParallelProjectionOn()

renwin = vtkRenderWindow()
renwin.AddRenderer(ren)

rwi = vtkRenderWindowInteractor()
rwi.SetRenderWindow(renwin)
#rwi.LightFollowCameraOn();
rwi.Initialize()
rwi.Start()


