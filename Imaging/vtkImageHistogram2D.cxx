#include "vtkImageHistogram2D.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageHistogram2D, "$Revision: 1.2 $");
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
    // store these ranges in our own ivars as well!
    double in1range[2], in2range[2];
    in1Data->GetScalarRange(in1range);
    in2Data->GetScalarRange(in2range);

    int input1bins = self->GetInput1Bins();
    int input2bins = self->GetInput2Bins();

    // clear output data
    memset(outData, 0, input1bins * input2bins);
        
    double in1binWidth = (in1range[1] - in1range[0]) / input1bins;
    double in2binWidth = (in2range[1] - in2range[0]) / input2bins;

    // calculate number of elements
    unsigned long noe = in1Data->GetDimensions()[0] *
        in1Data->GetDimensions()[1] * in1Data->GetDimensions()[2];

    T in1val, in2val;
    unsigned bin1, bin2;
    for (unsigned long i = 0; i < noe; i++)
    {
        in1val = *in1Ptr;
        in2val = *in2Ptr;
        
        bin1 = (unsigned)((in1val - in1range[0]) / in1binWidth);
        bin2 = (unsigned)((in2val - in2range[0]) / in2binWidth);

        // increment the correct bin
        *(outPtr + bin2 * input1bins + bin1) += 1;
    }
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

  // check origin, dimensions, spacing equal
  bool originEqual = true, spacingEqual = true, dimensionsEqual = true;
  for (int i = 0; i < 3 && originEqual && spacingEqual && dimensionsEqual; i++)
  {
      originEqual =
          this->GetInput1()->GetOrigin()[i] ==
          this->GetInput2()->GetOrigin()[i];

      spacingEqual =
          this->GetInput1()->GetSpacing()[i] ==
          this->GetInput2()->GetSpacing()[i];

      dimensionsEqual = 
          this->GetInput1()->GetDimensions()[i] ==
          this->GetInput2()->GetDimensions()[i];
      
  }

  if (! (originEqual && spacingEqual && dimensionsEqual))
  {
      vtkErrorMacro(<< "Execute: The two inputs have to match w.r.t. origin, "
                    "spacing and dimensions.");
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

