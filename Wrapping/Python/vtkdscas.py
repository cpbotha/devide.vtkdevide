"""
Code greatly simplified (more manual) so that McMillan installer has an
easier time finding the imports.

$Id: vtkdscas.py,v 1.1 2003/01/08 16:14:50 cpbotha Exp $
"""

import os

try:
    if os.name == 'posix':
        from libvtkdscasIOPython import *
        from libvtkdscasRenderingPython import *
    else:
        from vtkdscasIOPython import *
        from vtkdscasRenderingPython import *
except ImportError,e:
    print e
    pass
