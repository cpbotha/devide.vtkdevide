// vtkEventOutputWindow.cxx copyright (c) 2003 Charl P. Botha cpbotha@ieee.org
// $Id: vtkEventOutputWindow.cxx,v 1.1 2003/09/23 14:21:11 cpbotha Exp $
// vtkOutputWindow derivative that InvokeEvents instead of trying to display by itself

#include "vtkEventOutputWindow.h"
#include "vtkCommand.h"
#include "vtkObjectFactory.h"
#include "string.h"

vtkCxxRevisionMacro(vtkEventOutputWindow, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkEventOutputWindow);

#define MAX_TEXT_LEN 4095

vtkEventOutputWindow::vtkEventOutputWindow()
{
  this->TextType = 0;
  this->Text = new char[MAX_TEXT_LEN + 1];
  this->Text[MAX_TEXT_LEN] = (char)0;
}

vtkEventOutputWindow::~vtkEventOutputWindow()
{
  delete[] this->Text;
}

void vtkEventOutputWindow::SetupForEvent(const char *text)
{
  int tlen = strlen(text);
  // make sure we don't copy more than MAX_TEXT_LEN characters
  strncpy(this->Text, text, tlen < MAX_TEXT_LEN ? tlen : MAX_TEXT_LEN);
}

void vtkEventOutputWindow::DisplayText(const char *text)
{
  if (!text)
    {
    return;
    }

  this->SetupForEvent(text);
  // set us up for just text
  this->TextType = 0;
  // finally invoke the event
  this->InvokeEvent(vtkCommand::WarningEvent);
}

void vtkEventOutputWindow::DisplayErrorText(const char *text)
{
  if (!text)
    {
    return;
    }

  this->SetupForEvent(text);
  // set us up for just text
  this->TextType = 1;
  // finally invoke the event
  this->InvokeEvent(vtkCommand::ErrorEvent);
}

void vtkEventOutputWindow::DisplayWarningText(const char *text)
{
  if (!text)
    {
    return;
    }

  this->SetupForEvent(text);
  // set us up for just text
  this->TextType = 2;
  // finally invoke the event
  this->InvokeEvent(vtkCommand::WarningEvent);
}

void vtkEventOutputWindow::DisplayGenericWarningText(const char *text)
{
  if (!text)
    {
    return;
    }

  this->SetupForEvent(text);
  // set us up for just text
  this->TextType = 3;
  // finally invoke the event
  this->InvokeEvent(vtkCommand::WarningEvent);

}

void vtkEventOutputWindow::DisplayDebugText(const char *text)
{
  if (!text)
    {
    return;
    }

  this->SetupForEvent(text);
  // set us up for just text
  this->TextType = 4;
  // finally invoke the event
  this->InvokeEvent(vtkCommand::WarningEvent);

}


