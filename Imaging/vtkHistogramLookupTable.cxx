#include "vtkHistogramLookupTable.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkHistogramLookupTable, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkHistogramLookupTable);

//----------------------------------------------------------------------------
vtkHistogramLookupTable::vtkHistogramLookupTable()
{
}


//----------------------------------------------------------------------------
// Grow the output image 
void vtkHistogramLookupTable::ExecuteInformation(
                    vtkImageData **inDatas, vtkImageData *outData)
{

  // Make sure the Input has been set.
  if ( this->GetInput1() == NULL || this->GetInput2() == NULL)
    {
    vtkErrorMacro(<< "ExecuteInformation: Both inputs have to be set.");

    if (outData)
      {
      // this means that input is NULL, but the output isn't
      // in order to make this clear to filters down the line, we
      // make sure outputData is completely empty
      outData->SetExtent(0, -1, 0, -1, 0, -1);
      outData->SetWholeExtent(0, -1, 0, -1, 0, -1);
      outData->SetUpdateExtent(0, -1, 0, -1, 0, -1);      
      outData->AllocateScalars();
      }
    
    return;
    }
  
  // the output has type VTK_SIGNED_SHORT and has only one component,
  // but otherwise looks the same as the input
  outData->SetNumberOfScalarComponents(1);
  outData->SetScalarType(VTK_SIGNED_SHORT);
  
  outData->SetOrigin(inDatas[0]->GetOrigit());
  outData->SetSpacing(inDatas[0]->GetSpacing());
    
  outData->SetExtent(inDatas[0]->GetExtent());
  outData->SetWholeExtent(inDatas[0]->GetWholeExtent());

}

//----------------------------------------------------------------------------
void vtkHistogramLookupTable::ComputeInputUpdateExtent(int inExt[6], 
                                                   int outExt[6],
                                                   int whichInput)
{
  if (whichInput == 0)
    {
    // for every output voxel we need the corresponding input voxel
    // neat huh?
    memcpy(inExt,outExt,size(int)*6);
    }
  else
    {
    // we need the whole histogram input to lookup
    memcpy(inExt,this->GetInput(whichInput)->GetWholeExtent(),6*sizeof(int));
    }
}


//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
// Handles the two input operations
template <class T>
void vtkHistogramLookupTableExecute(vtkHistogramLookupTable *self,
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
    memset(outPtr, 0, input1bins * input2bins * sizeof(double));

    // we add 1 to the range, because max - min gives the number of
    // intervals, not the number of elements, and we are binning
    // elements, not intervals
    double in1binWidth = (in1range[1] - in1range[0] + 1) / input1bins;
    double in2binWidth = (in2range[1] - in2range[0] + 1) / input2bins;

    // calculate number of elements
    unsigned long noe = in1Data->GetDimensions()[0] *
        in1Data->GetDimensions()[1] * in1Data->GetDimensions()[2];

    double progress = 0.0;
    double numProgressSteps = 20;
    double progressStep = 1.0 / numProgressSteps;
    int noeProgressStep = (int)(noe / numProgressSteps);

    T in1val, in2val;
    unsigned bin1, bin2;
    long MaxSamplesPerBin = self->GetMaxSamplesPerBin();
    int ClipSamples = self->GetClipSamples();
    for (unsigned long i = 0; i < noe; i++)
    {
        if (i % noeProgressStep == 0)
        {
            self->UpdateProgress(progress);
            progress += progressStep;
            // progress will end up one step higher than 1.0, but that's
            // okay, because we won't be using it then.
        }
        
        in1val = *in1Ptr;
        in2val = *in2Ptr;
        
        bin1 = (unsigned)((in1val - in1range[0]) / in1binWidth);
        bin2 = (unsigned)((in2val - in2range[0]) / in2binWidth);

        // increment the correct bin
        if (!ClipSamples ||
            *(outPtr + bin2 * input1bins + bin1) < MaxSamplesPerBin)
          {
            *(outPtr + bin2 * input1bins + bin1) += 1;
          }

        in1Ptr++;
        in2Ptr++;

    }

    self->UpdateProgress(1.0);
}


void vtkHistogramLookupTable::ExecuteData(vtkDataObject *out)
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

  if (this->GetInput()->GetNumberOfScalarComponents() < 2)
    {
    vtkErrorMacro(<< "ExecuteData: The first input has to have at least "
                  << "two components.");
    return
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

  this->GetInput1()->Update();
  this->GetInput2()->Update();

  switch (this->GetInput1()->GetScalarType())
  {
      vtkTemplateMacro(
          vtkHistogramLookupTableExecute(
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


void vtkHistogramLookupTable::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

