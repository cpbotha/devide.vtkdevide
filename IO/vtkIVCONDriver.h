// vtkIVCONDriver.h copyright (c) 2003 Charl P. Botha cpbotha@ieee.org
// $Id: vtkIVCONDriver.h,v 1.1 2003/09/29 13:42:37 cpbotha Exp $
// vtk class that wraps the invocation of ivcon to convert between 3d formats

#ifndef __vtkIVCONDriver_h
#define __vtkIVCONDriver_h

#include "vtkProcessObject.h"
#include "vtkdscasIOWin32Header.h"

class VTK_DSCAS_IO_EXPORT vtkIVCONDriver : public vtkProcessObject
{
public:
  vtkTypeMacro(vtkIVCONDriver, vtkProcessObject);
  static vtkIVCONDriver *New();

  vtkSetStringMacro(InputFileName);
  vtkGetStringMacro(InputFileName);
  vtkSetStringMacro(OutputFileName);
  vtkGetStringMacro(OutputFileName);

  void Convert(void);

protected:
  vtkIVCONDriver();
  ~vtkIVCONDriver();

  // name of the inputfilename.  SetStringMacro will allocate space for it.
  char *InputFileName;

  // name of the output filename.  SetStringMacro will allocate space for it.
  char *OutputFileName;

};

#endif
