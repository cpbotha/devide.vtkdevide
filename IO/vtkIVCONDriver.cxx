#include "vtkIVCONDriver.h"
#include "vtkObjectFactory.h"
#include "ivcon/ivconv.h"

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
  if (this->InputFileName && strlen(this->InputFileName) > 0 && 
      this->OutputFileName && strlen(this->OutputFileName) > 0)
    {
    if (this->GetMTime() > this->ConvertTime.GetMTime())
      {
        // do the actual conversion
        IVCONV ivconv;
        bool ivconvOk;
        ivconvOk = ivconv.read(this->InputFileName);
        if (!ivconvOk)
          {
          vtkErrorMacro(<<"Converter code could not read input file " << this->InputFileName);
          return;
          }

        ivconvOk = ivconv.write(this->OutputFileName);
        if (!ivconvOk)
          {
          vtkErrorMacro(<<"Converter code could not write output file " << this->OutputFileName);
          return;
          }

        // update timestamp
        this->ConvertTime.Modified();
      }
    }
  else
    {
    vtkErrorMacro(<< "Please set InputFileName and OutputFileName before calling Convert.");
    }
}

