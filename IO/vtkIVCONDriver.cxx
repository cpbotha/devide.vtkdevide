#include "vtkIVCONDriver.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkIVCONDriver);

vtkIVCONDriver::vtkIVCONDriver() : vtkProcessObject()
{
  this->InputFileName = NULL;
  this->OutputFileName = NULL;
}

vtkIVCONDriver::~vtkIVCONDriver()
{
  if (this->InputFileName)
    {
    delete[] this->InputFileName;
    }
  if (this->OutputFileName)
    {
    delete[] this->OutputFileName;
    }
}

void vtkIVCONDriver::Convert(void)
{
}

