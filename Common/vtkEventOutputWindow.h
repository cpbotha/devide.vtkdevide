// vtkEventOutputWindow.h copyright (c) 2003 Charl P. Botha cpbotha@ieee.org
// $Id: vtkEventOutputWindow.h,v 1.1 2003/09/23 14:21:11 cpbotha Exp $
// vtkOutputWindow derivative that InvokeEvents instead of trying to display by itself

#ifndef __vtkEventOutputWindow_h
#define __vtkEventOutputWindow_h

// .NAME vtkEventOutputWindow - output is resent as an event.
// .SECTION Description
// Sends all output to any registered observers as events.  Errors are sent as
// ErrorEvents and everything else as WarningEvents.  This instance can be queried
// for the actual output text (and maybe later for more information).  The TextType
// can be used to distinguish between the different output types.
// To use this class, instantiate it and then call SetInstance(this).
// 


#include "vtkOutputWindow.h"
#include "vtkdscasCommonWin32Header.h"

class VTK_DSCAS_COMMON_EXPORT vtkEventOutputWindow : public vtkOutputWindow
{
public:
  static vtkEventOutputWindow *New();
  vtkTypeRevisionMacro(vtkEventOutputWindow,vtkOutputWindow);

  virtual void DisplayText(const char*);
  virtual void DisplayErrorText(const char*);
  virtual void DisplayWarningText(const char*);
  virtual void DisplayGenericWarningText(const char*);
  virtual void DisplayDebugText(const char*);

  const char *GetText(void)
  {
    return this->Text;
  }

  // Description: TextType indicates the type of output.
  // 0: Text, 1: ErrorText, 2: WarningText, 3: GenericWarningText, 4: DebugText
  vtkGetMacro(TextType, int);

protected:
  vtkEventOutputWindow();
  ~vtkEventOutputWindow();
  void SetupForEvent(const char *);
  char *Text;
  int TextType;

private:
  vtkEventOutputWindow(const vtkEventOutputWindow&);  //Not implemented
  void operator=(const vtkEventOutputWindow&);  //Not implemented
};

#endif
