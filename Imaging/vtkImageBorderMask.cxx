#include "vtkImageBorderMask.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageBorderMask, "$Revision: 1.4 $");
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
                               int Borders[],
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

  // now do sanity checking on Borders
  for (int i = 0; i < 6; i++)
    {
    if (Borders[i] < 0)
      {
      vtkGenericWarningMacro(<< "Execute: " << i << "'th border negative.  "
                             << "Assuming default 1.");
      Borders[i] = 1;
      }
    }

  double scalarRange[2];
  input->GetScalarRange(scalarRange);
  IT inputMinimum = (IT)(scalarRange[0]);
  IT inputMaximum = (IT)(scalarRange[1]);
  
  int zOffset, yOffset, xOffset;
  bool InZBorder, InYBorder, InXBorder;

  for (int z = 0; z < dims[2]; z++)
    {
    zOffset = z * dims[1] * dims[0];

    if (z < Borders[4] || z >= dims[2] - Borders[5])
      {
      InZBorder = true;
      }
    else
      {
      InZBorder = false;
      }

    // now the interior y-lines (plus border x left and right)
    for (int y = 0; y < dims[1]; y++)
      {
      yOffset = zOffset + y * dims[0];

      if (y < Borders[2] || y >= dims[1] - Borders[3])
        {
        InYBorder = true;
        }
      else
        {
        InYBorder = false;
        }

      for (int x = 0; x < dims[0]; x++)
        {
        xOffset = yOffset + x;

        if (x < Borders[0] || x >= dims[0] - Borders[1])
          {
          InXBorder = true;
          }
        else
          {
          InXBorder = false;
          }

        if (InZBorder || InYBorder || InXBorder)
          {
          outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
            BorderMode, BorderValue,
            inPtr[xOffset], inputMinimum, inputMaximum);
          }
        else
          {
          outPtr[xOffset] = vtkImageBorderMaskReplaceValue(
            InteriorMode, InteriorValue,
            inPtr[xOffset], inputMinimum, inputMaximum);
          }
          
        } // for (int x ...

      } // for (int y = 0 ...
    } // for (int z = 0 ...
}

vtkImageBorderMask::vtkImageBorderMask() : vtkSimpleImageToImageFilter()
{
  // use BorderValue for border
  this->SetBorderMode(0);
  this->SetBorderValue(1);
  this->SetBorders(1,1,1,1,1,1);

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
    vtkTemplateMacro9(vtkImageBorderMaskExecute, input, output,
                      (VTK_TT *)(inPtr), (VTK_TT *)(outPtr),
                      this->BorderMode, this->BorderValue,
                      this->Borders,
                      this->InteriorMode, this->InteriorValue);
    default:
      vtkGenericWarningMacro("Execute: Unknown input ScalarType");
      return;
    }
}
