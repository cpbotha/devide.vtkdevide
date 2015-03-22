// vtkEventOutputWindow.cxx copyright (c) 2003 Charl P. Botha cpbotha@ieee.org
// $Id$
// vtkOutputWindow derivative that InvokeEvents instead of trying to display by itself

#include "vtkEventOutputWindow.h"
#include "vtkCommand.h"
#include "vtkObjectFactory.h"
#include "string.h"

vtkCxxRevisionMacro(vtkEventOutputWindow, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkEventOutputWindow);

#define MAX_TEXT_LEN 4095

vtkEventOutputWindow::vtkEventOutputWindow()
{
  this->TextType = 0;
  this->Text = new char[MAX_TEXT_LEN + 1];
  this->Text[MAX_TEXT_LEN] = this->Text[0] = (char)0;
}

vtkEventOutputWindow::~vtkEventOutputWindow()
{
  delete[] this->Text;
}

void vtkEventOutputWindow::SetupForEvent(const char *text)
{
  int tlen = strlen(text) + 1; // add one for null termination
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
  this->InvokeEvent(vtkCommand::WarningEvent, (void*)text);
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
  this->InvokeEvent(vtkCommand::ErrorEvent, (void*)text);
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
  this->InvokeEvent(vtkCommand::WarningEvent, (void*)text);
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
  this->InvokeEvent(vtkCommand::WarningEvent, (void*)text);

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
  this->InvokeEvent(vtkCommand::WarningEvent, (void*)text);

}


