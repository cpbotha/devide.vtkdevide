"""
Code greatly simplified (more manual) so that McMillan installer has an
easier time finding the imports.

$Id: vtkdevide.py,v 1.3 2004/11/01 11:09:05 cpbotha Exp $
"""

import os

try:
    if os.name == 'posix':
        from libvtkdevideCommonPython import *        
        from libvtkdevideHybridPython import *
        from libvtkdevideImagingPython import *
        from libvtkdevideIOPython import *
        from libvtkdevideRenderingPython import *
        from libvtkdevideExternalPython import *
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
