#include "vtkImageGreyscaleReconstruct3D.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageGreyscaleReconstruct3D, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImageGreyscaleReconstruct3D);

//----------------------------------------------------------------------------
vtkImageGreyscaleReconstruct3D::vtkImageGreyscaleReconstruct3D()
{
  this->SetDual(0);
}


//----------------------------------------------------------------------------
// Grow the output image 
void vtkImageGreyscaleReconstruct3D::ExecuteInformation(
                    vtkImageData **inDatas, vtkImageData *outData)
{
    inDatas[1]->UpdateInformation();
    // copy all metadata (I hope) from the J marker image
    outData->CopyStructure(inDatas[1]);
}

//----------------------------------------------------------------------------
void vtkImageGreyscaleReconstruct3D::ComputeInputUpdateExtent(int inExt[6], 
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
void vtkImageGreyscaleReconstruct3DExecute(vtkImageGreyscaleReconstruct3D *self,
                                vtkImageData *in1Data, T *in1Ptr,
                                vtkImageData *in2Data, T *in2Ptr,
                                vtkImageData *outData)
{
  int Dual = self->GetDual();

  // first thing, copy J (second input) over as output (we'll work in
  // the output)
  outData->DeepCopy(in2Data);
  int xdim = outData->GetDimensions()[0];
  int ydim = outData->GetDimensions()[1];
  int zdim = outData->GetDimensions()[2];

  // Nplus neighbourhood pointers will be stored in this
  T *NpPtrs[13];
  // these are the Npi offsets we need to extract the neighbourhood
  int Npi[13];

  Npi[0] = -ydim*xdim - xdim - 1; // x-1, y-1, z-1
  Npi[1] = -ydim*xdim - xdim;     // x,   y-1, z-1
  Npi[2] = -ydim*xdim - xdim +1;  // x+1, y-1, z-1

  Npi[3] = -ydim*xdim - 1;        // x-1, y,   z-1
  Npi[4] = -ydim*xdim;            // x,   y,   z-1
  Npi[5] = -ydim*xdim + 1;        // x+1, y,   z-1

  Npi[6] = -ydim*xdim + xdim - 1; // x-1, y+1, z-1
  Npi[7] = -ydim*xdim + xdim;     // x,   y+1, z-1
  Npi[8] = -ydim*xdim + xdim + 1; // x+1, y+1, z-1

  Npi[9] = -xdim - 1;             // x-1, y-1, z
  Npi[10] = -xdim;                // x,   y-1, z
  Npi[11] = -xdim + 1;            // x+1, y-1, z

  Npi[12] = -1;                   // x-1, y,   z

  T *pPtr;
  T *Iptr;
  T minNp;
  int i;
  if (Dual)
    {
    // 1. scan D1 in raster order
    pPtr = (T*)(outData->GetScalarPointer());
    Iptr = in1Ptr;
    for (int z = 0; z < zdim; z++)
      {
      for (int y = 0; y < ydim; y++)
        {
        for (int x = 0; x < xdim; x++)
          {
          // J(p) max(min(J(q), q E N+(p) U p) ,I(p))
          // i.e. find the minimum in the N+ neighbourhood of P and P
          // replace J(p = x,y,z) with the max of that minimum and
          // I(p)

          // a. determine minimum of Np neighbourhood
          for (i = 0; i < 13; i++)
            {
            // calculate modified pointer
            NpPtrs[i] = pPtr + Npi[i];
            }

          // now zero things which are over the boundary
          if (x == 0)
            {
            NpPtrs[0] = NpPtrs[3] = NpPtrs[6] = NpPtrs[9] = NpPtrs[12] = \
              (T*)NULL;
            }
          else if (x == xdim - 1)
            {
            NpPtrs[2] = NpPtrs[5] = NpPtrs[8] = NpPtrs[11] = (T*)NULL;
            }
          if (y == 0)
            {
            NpPtrs[0] = NpPtrs[1] = NpPtrs[2] = NpPtrs[9] = NpPtrs[10] = \
              NpPtrs[11] = (T*)NULL;
            }
          else if (y == ydim - 1)
            {
            NpPtrs[6] = NpPtrs[7] = NpPtrs[8] = (T*)NULL;
            }
          if (z == 0)
            {
            NpPtrs[0] = NpPtrs[1] = NpPtrs[2] = \
              NpPtrs[3] = NpPtrs[4] = NpPtrs[5] =                \
              NpPtrs[6] = NpPtrs[7] = NpPtrs[8] = (T*)NULL;
            }
          // there are no z+1 cases

          // now determine the minimum
          minNp = *pPtr;
          for (i = 0; i < 13; i++)
            {
            if (NpPtrs[i] && *NpPtrs[i] < minNp)
              {
              minNp = *NpPtrs[i];
              }
            }
          *pPtr = minNp > *Iptr ? minNp : *Iptr;
          
          // increment pointer
          pPtr++;
          Iptr++;
          }
        }
      }
    // 2. scan D1 in anti-raster order

    // 3. propagation step
    
    } // if (Dual) ...
  else
    {
    // not implemented yet...
    }

    self->UpdateProgress(1.0);
}


void vtkImageGreyscaleReconstruct3D::ExecuteData(vtkDataObject *out)
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

  this->GetInput1()->Update();
  this->GetInput2()->Update();

  switch (this->GetInput2()->GetScalarType())
  {
      vtkTemplateMacro(
          vtkImageGreyscaleReconstruct3DExecute(
              this,
              this->GetInput1(), static_cast<VTK_TT*>(in1Ptr),
              this->GetInput2(), static_cast<VTK_TT*>(in2Ptr),
              output
              )
          );

    default:                              
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }


  
}


void vtkImageGreyscaleReconstruct3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Dual: " << this->Dual << "\n";
}



