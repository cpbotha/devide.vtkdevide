"""
Code greatly simplified (more manual) so that McMillan installer has an
easier time finding the imports.

$Id: vtkdevide.py,v 1.1 2004/01/15 11:00:57 cpbotha Exp $
"""

import os

try:
    if os.name == 'posix':
        from libvtkdevideCommonPython import *        
        from libvtkdevideHybridPython import *        
        from libvtkdevideIOPython import *
        from libvtkdevideRenderingPython import *
    else:
        from vtkdevideCommonPython import *        
        from vtkdevideHybridPython import *        
        from vtkdevideIOPython import *
        from vtkdevideRenderingPython import *
except ImportError,e:
    print e
    pass
