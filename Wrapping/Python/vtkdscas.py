"""
Code greatly simplified (more manual) so that McMillan installer has an
easier time finding the imports.

$Id: vtkdscas.py,v 1.3 2003/09/23 14:35:52 cpbotha Exp $
"""

import os

try:
    if os.name == 'posix':
        from libvtkdscasCommonPython import *        
        from libvtkdscasHybridPython import *        
        from libvtkdscasIOPython import *
        from libvtkdscasRenderingPython import *
    else:
        from vtkdscasCommonPython import *        
        from vtkdscasHybridPython import *        
        from vtkdscasIOPython import *
        from vtkdscasRenderingPython import *
except ImportError,e:
    print e
    pass
