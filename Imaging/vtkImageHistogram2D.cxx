#include "vtkImageHistogram2D.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageHistogram2D, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImageHistogram2D);

//----------------------------------------------------------------------------
vtkImageHistogram2D::vtkImageHistogram2D()
{
    this->SetInput1Bins(256);
    this->SetInput2Bins(256);
}


//----------------------------------------------------------------------------
// Grow the output image 
void vtkImageHistogram2D::ExecuteInformation(
                    vtkImageData **vtkNotUsed(inDatas), vtkImageData *outData)
{
    outData->SetNumberOfScalarComponents(1);
    outData->SetScalarType(VTK_DOUBLE);
    outData->SetOrigin(0,0,0);
    outData->SetSpacing(1,1,1);
    outData->SetExtent(0, this->Input1Bins - 1, 0, this->Input2Bins - 1, 0, 1);
    outData->SetWholeExtent(0, this->Input1Bins - 1,
                            0, this->Input2Bins - 1, 0, 1);
}

//----------------------------------------------------------------------------
void vtkImageHistogram2D::ComputeInputUpdateExtent(int inExt[6], 
                                                   int outExt[6],
                                                   int whichInput)
{
    // get the whole image for input 2
    memcpy(inExt,this->GetInput(whichInput)->GetWholeExtent(),6*sizeof(int));
}


//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
// Handles the two input operations
template <class T>
void vtkImageHistogram2DExecute(vtkImageHistogram2D *self,
                                vtkImageData *in1Data, T *in1Ptr,
                                vtkImageData *in2Data, T *in2Ptr,
                                vtkImageData *outData, double *outPtr)
{
}


void vtkImageHistogram2D::ExecuteData(vtkDataObject *out)
{
  // Make sure the Input has been set.
  if ( this->GetInput1() == NULL || this->GetInput2() == NULL)
    {
    vtkErrorMacro(<< "ExecuteData: Both inputs have to be set.");
    return;
    }
  
  // Too many filters have floating point exceptions to execute
  // with empty input/ no request.
  if (this->UpdateExtentIsEmpty(out))
    {
    return;
    }

  if (this->GetInput1()->GetScalarType() != this->GetInput2()->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, " <<
    this->GetInput1()->GetScalarType() << " and input2 ScalarType " <<
    this->GetInput2()->GetScalarType() << ", should match");
    return;
    }

  // get metadata across
  this->ExecuteInformation();
  // make sure output has been allocated
  vtkImageData *output = vtkImageData::SafeDownCast(out);
  output->AllocateScalars();

  
  // now let's start with the actual work
  void *in1Ptr = this->GetInput1()->GetScalarPointer();
  void *in2Ptr = this->GetInput2()->GetScalarPointer();
  double *outPtr = (double*)output->GetScalarPointer();

  switch (this->GetInput1()->GetScalarType())
  {
      vtkTemplateMacro(
          vtkImageHistogram2DExecute(
              this,
              this->GetInput1(), static_cast<VTK_TT*>(in1Ptr),
              this->GetInput2(), static_cast<VTK_TT*>(in2Ptr),
              output, outPtr
              )
          );

    default:                              
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }


  
}


void vtkImageHistogram2D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Input1Bins: " << this->Input1Bins << "\n";
  os << indent << "Input2Bins: " << this->Input2Bins << "\n";  
  
}

