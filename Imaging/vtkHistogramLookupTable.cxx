#include "vtkHistogramLookupTable.h"

#include "vtkImageData.h"
#include "vtkImageProgressIterator.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkHistogramLookupTable, "$Revision: 1.5 $");
vtkStandardNewMacro(vtkHistogramLookupTable);

#define OUTPUT_TYPE signed short

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
  
  // the output has type VTK_SHORT and has only one component,
  // but otherwise looks the same as the input
  outData->SetNumberOfScalarComponents(1);
  outData->SetScalarType(VTK_SHORT);
  
  outData->SetOrigin(inDatas[0]->GetOrigin());
  outData->SetSpacing(inDatas[0]->GetSpacing());
    
  outData->SetExtent(outData->GetUpdateExtent());
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
    memcpy(inExt,outExt,sizeof(int)*6);
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
template <class I1T, class I2T>
void vtkHistogramLookupTableExecute(vtkHistogramLookupTable *self,
                                    vtkImageData *in1Data,
                                    vtkImageData *in2Data,
                                    vtkImageData *outData, int outExt[6],
                                    I1T *, I2T *)

{
  in1Data->Update();
  in2Data->Update();
  
  vtkImageIterator<I1T> in1It(in1Data, outExt);

  // we're single-threaded, so we can pass thread-id 0
  vtkImageProgressIterator<OUTPUT_TYPE> outIt(outData, outExt, self, 0);
  I1T tempVal1, tempVal2;

  int histogramExtent[6];
  in2Data->GetWholeExtent(histogramExtent);
  double histogramSpacing[3];
  in2Data->GetSpacing(histogramSpacing);
  double histogramOrigin[3];
  in2Data->GetOrigin(histogramOrigin);
  
  int bins1 = histogramExtent[1] - histogramExtent[0] + 1;
  int bins2 = histogramExtent[3] - histogramExtent[2] + 1;
  double binWidth1 = histogramSpacing[0];
  double binWidth2 = histogramSpacing[1];
  double range1min = histogramOrigin[0];
  double range2min = histogramOrigin[1];

  unsigned int bin1, bin2;
  I2T *in2Ptr = (I2T *)(in2Data->GetScalarPointer());
  
  while (!outIt.IsAtEnd())
    {
    I1T *inSI = in1It.BeginSpan();
    OUTPUT_TYPE *outSI = outIt.BeginSpan();
    OUTPUT_TYPE *outSIEnd = outIt.EndSpan();
    while (outSI < outSIEnd)
      {
      tempVal1 = *inSI;
      ++inSI; // go to next scalar component
      tempVal2 = *inSI;
      ++inSI; // go to next pair

      // now check where in the histogram tempVal1 and tempVal2 are
      bin1 = (unsigned)((tempVal1 - range1min) / binWidth1);
      bin2 = (unsigned)((tempVal2 - range2min) / binWidth2);
      if (*(in2Ptr + bin2 * bins1 + bin1) > 0)
        {
        *outSI = 1;
        }
      else
        {
        *outSI = 0;
        }
      ++outSI;
      } // while (outSI < outSIEnd)
    in1It.NextSpan();
    outIt.NextSpan();
    } // while (!outIt.IsAtEnd()) ...

    self->UpdateProgress(1.0);
}

template <class T>
void vtkHistogramLookupTableExecute1(vtkHistogramLookupTable *self,
                                     vtkImageData *input1,
                                     vtkImageData *input2,
                                     vtkImageData *output, int outExt[6],
                                     T *)
{
  // second stage of type specialization
  switch (input2->GetScalarType())
    {
    vtkTemplateMacro7(vtkHistogramLookupTableExecute, self,
                      input1, input2, output, outExt,
                      static_cast<T *>(0), static_cast<VTK_TT *>(0));
    default:
      vtkGenericWarningMacro(<< "Execute: Unknown ScalarType for second input.");
      return;
    }
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

  this->GetInput1()->Update();
  this->GetInput2()->Update();

  int outExt[6];
  output->GetUpdateExtent(outExt);
  
  switch (this->GetInput1()->GetScalarType())
    {
    // first part of a two step specialization - first we do the
    // input1 type
    vtkTemplateMacro6(vtkHistogramLookupTableExecute1, this,
                      this->GetInput1(), this->GetInput2(),
                      output, outExt,
                      static_cast<VTK_TT *>(0));

    default:                              
      vtkErrorMacro(<< "Execute: Unknown ScalarType for first input.");
      return;
    }


  
}


void vtkHistogramLookupTable::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

