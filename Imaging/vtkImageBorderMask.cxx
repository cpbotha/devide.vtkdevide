#include "vtkImageBorderMask.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageBorderMask, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkImageBorderMask);



template <class IT>
inline IT vtkImageBorderMaskReplaceValue(int mode, int defaultvalue,
                                         IT inputvalue,
                                         IT inputMinimum, IT inputMaximum)
{
  IT ReplaceValue;
  
  if (mode == 0)
    {
    ReplaceValue = (IT)defaultvalue;
    }
  else if (mode == 1)
    {
    ReplaceValue = inputvalue;
    }
  else if (mode == 2)
    {
    ReplaceValue = inputMinimum;
    }
  else
    {
    ReplaceValue = inputMaximum;
    }

  return ReplaceValue;
}


// The switch statement in Execute will call this method with
// the appropriate input type (IT). Note that this example assumes
// that the output data type is the same as the input data type.
// This is not always the case.
template <class IT>
void vtkImageBorderMaskExecute(vtkImageData* input,
                               vtkImageData* output,
                               IT* inPtr, IT* outPtr,
                               int BorderMode, int BorderValue,
                               int InteriorMode, int InteriorValue)
{
  if (input->GetScalarType() != output->GetScalarType())
    {
    vtkGenericWarningMacro(<< "Execute: input ScalarType, "
                           << input->GetScalarType()
                           << ", must match out ScalarType "
                           << output->GetScalarType());
    return;
    }

  // get the input up to date
  input->Update();
  

  int dims[3];
  input->GetDimensions(dims);

  // guarantee at least 3x3x1 volume
  if (dims[0] < 3 || dims[1] < 3)
    {
    vtkGenericWarningMacro(<< "Execute: Input image/volume too small to "
                           << "create border mask.");
    return;
    }

  double scalarRange[2];
  input->GetScalarRange(scalarRange);
  IT inputMinimum = (IT)(scalarRange[0]);
  IT inputMaximum = (IT)(scalarRange[1]);
  
  int zintStart, zintEnd;
  if (dims[2] < 3)
    {
    // we're either looking at an image or a 2-slice volume
    zintStart = 0;
    zintEnd = dims[2];
    }
  else
    {
    // we have a real volume
    zintStart = 1;
    zintEnd = dims[2] - 1;
    }

  int zOffset, yOffset, xOffset;
  int topAndBottomY[2];

  for (int z = zintStart; z < zintEnd; z++)
    {
    zOffset = z * dims[1] * dims[0];

	// do top and bottom y-rows
	topAndBottomY[0] = 0;
	topAndBottomY[1] = dims[1] - 1;
	for (int i = 0; i < 2; i++)
	  {
	  yOffset = zOffset + topAndBottomY[i] * dims[0];

	  for (int x = 0; x < dims[0]; x++)
        {
        xOffset = yOffset + x;
        outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
              BorderMode, BorderValue,
              inPtr[xOffset], inputMinimum, inputMaximum);
          
        } // for (int x = 0 ...
      } // for (int i = 0 ...

    // now the interior y-lines (plus border x left and right)
    for (int y = 1; y < dims[1] - 1; y++)
      {
      yOffset = zOffset + y * dims[0];

      // leftmost x BOUNDARY
      xOffset = yOffset + 0;
      outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
	    		BorderMode, BorderValue,
		    	inPtr[xOffset], inputMinimum, inputMaximum);
      
      // INTERIOR
      for (int x = 1; x < dims[0] - 1; x++)
        {
        xOffset = yOffset + x;

        outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
          InteriorMode, InteriorValue,
          inPtr[xOffset], inputMinimum, inputMaximum);
          
        }

      // rightmost x BOUNDARY
      xOffset = yOffset + dims[0] - 1;
      outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
        BorderMode, BorderValue,
        inPtr[xOffset], inputMinimum, inputMaximum);
      
      } // for (int y = 1 ...
    } // for (int z = zintStart ...

  // then let's do the z==0 and z==dims[2] - 1 boundaries (if
  // applicable)
  if (zintStart == 1)
    {
    int BottomAndTopZ[2];
    BottomAndTopZ[0] = 0;
    BottomAndTopZ[1] = dims[2] - 1;

    for (int i = 0; i < 2; i++)
      {
      zOffset = BottomAndTopZ[i] * dims[1] * dims[0];
      
      for (int y = 0; y < dims[1]; y++)
        {
        yOffset = zOffset + y * dims[0];
        for (int x = 0; x < dims[0]; x++)
          {
          xOffset = yOffset + x;
          
          outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
            BorderMode, BorderValue,
            inPtr[xOffset], inputMinimum, inputMaximum);
          }
        }
      } // for (int i = 0 ...
    } // if (zintStart == 1 ...
}

vtkImageBorderMask::vtkImageBorderMask() : vtkSimpleImageToImageFilter()
{
  // use BorderValue for border
  this->SetBorderMode(0);
  this->SetBorderValue(1);

  // user InteriorValue for interior
  this->SetInteriorMode(0);
  this->SetInteriorValue(0);
}

void vtkImageBorderMask::SimpleExecute(vtkImageData* input,
                                       vtkImageData* output)
{

  void* inPtr = input->GetScalarPointer();
  void* outPtr = output->GetScalarPointer();

  switch(output->GetScalarType())
    {
    // This is simple a #define for a big case list. It handles
    // all data types vtk can handle.
    vtkTemplateMacro8(vtkImageBorderMaskExecute, input, output,
                      (VTK_TT *)(inPtr), (VTK_TT *)(outPtr),
                      this->BorderMode, this->BorderValue,
                      this->InteriorMode, this->InteriorValue);
    default:
      vtkGenericWarningMacro("Execute: Unknown input ScalarType");
      return;
    }
}
