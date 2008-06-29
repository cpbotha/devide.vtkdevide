"""
Module to import vtkdevide libraries.
"""

import os

# change this file and commit if you want this updated
stamp = '20080629-2155'
version = "$Revision$"

try:
    if os.name == 'posix':

        # took DL handling from Mathieu Malaterre's vtkgdcm.py

        # extremely important !
        # http://gcc.gnu.org/faq.html#dso
        # http://mail.python.org/pipermail/python-dev/2002-May/023923.html
        # http://wiki.python.org/moin/boost.python/CrossExtensionModuleDependencies
        # This is now merged in VTK 5.2:
        # http://vtk.org/cgi-bin/viewcvs.cgi/Wrapping/Python/vtk/__init__.py?r1=1.13&r2=1.14
        import sys
        orig_dlopen_flags = sys.getdlopenflags()
        try:
          import dl
        except ImportError:
          # are we on AMD64 ?
          try:
            import DLFCN as dl
          except ImportError:
            print "Could not import dl"
            dl = None
        if dl:
          #print "dl was imported"
          #sys.setdlopenflags(dl.RTLD_LAZY|dl.RTLD_GLOBAL)
          sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
  
        from libvtkdevideCommonPython import *        
        from libvtkdevideHybridPython import *
        from libvtkdevideImagingPython import *
        from libvtkdevideIOPython import *
        from libvtkdevideRenderingPython import *
        from libvtkdevideExternalPython import *

        # revert:
        sys.setdlopenflags(orig_dlopen_flags)
        del sys, dl
        del orig_dlopen_flags
  
    else:
        from vtkdevideCommonPython import *        
        from vtkdevideHybridPython import *
        from vtkdevideImagingPython import *
        from vtkdevideIOPython import *
        from vtkdevideRenderingPython import *
        from vtkdevideExternalPython import *
except ImportError,e:
    print e
    pass
