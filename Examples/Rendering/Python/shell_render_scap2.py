# $Id: shell_render_scap2.py,v 1.2 2004/01/15 11:00:53 cpbotha Exp $
# example to test shell renderer (*shudder*)

from vtk import *
from vtkdevide import *
import os
import stat

def ce_cb(obj, evt_name):
    if obj.GetKeyCode() == 'm':
        crm = splatmapper.GetRenderMode()
	crm = crm + 1
	if crm > 2:
	    crm = 0
        splatmapper.SetRenderMode(crm)
        print "rendermode switched to %d" % (crm)

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



dicom_dir = 'scap2'
filenames_init = os.listdir(dicom_dir)
filenames_init.sort() # make it alphabetical
# go through list of files in directory, perform trivial tests
# and create a new list of files 
dicom_fullnames = []

for filename in filenames_init:
    # make full filename
    fullname = os.path.join(dicom_dir, filename)
    # at the moment, we check that it's a regular file
    if stat.S_ISREG(os.stat(fullname)[stat.ST_MODE]):
	dicom_fullnames.append(fullname)
	
dicomr = vtkDICOMVolumeReader()
	
for fullname in dicom_fullnames:
    # this will simply add a file to the buffer list of the
    # vtkDICOMVolumeReader (will not set mtime)
    print "%s\n" % fullname
    dicomr.add_dicom_filename(fullname)


dicomr.SetSeriesInstanceIdx(0)

otf = vtkPiecewiseFunction()
otf.AddPoint(0.0, 0.0)
otf.AddPoint(99, 0.0)
otf.AddPoint(100, 1.0)
otf.AddPoint(2800.0, 1.0)

ctf = vtkColorTransferFunction()
ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(99, 0.0, 0.0, 0.0)
ctf.AddRGBPoint(100, 1.0, 0.937, 0.859)
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
splatmapper.SetInput(dicomr.GetOutput())
splatmapper.SetRenderMode(1)

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

