"""
Code greatly simplified (more manual) so that McMillan installer has an
easier time finding the imports.

$Id: vtkdscas.py,v 1.2 2003/06/30 09:59:37 cpbotha Exp $
"""

import os

try:
    if os.name == 'posix':
        from libvtkdscasHybridPython import *        
        from libvtkdscasIOPython import *
        from libvtkdscasRenderingPython import *
    else:
        from vtkdscasHybridPython import *        
        from vtkdscasIOPython import *
        from vtkdscasRenderingPython import *
except ImportError,e:
    print e
    pass
