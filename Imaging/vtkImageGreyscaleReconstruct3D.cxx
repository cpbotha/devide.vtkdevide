#include "vtkImageGreyscaleReconstruct3D.h"

#include <vtkstd/queue>
#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageGreyscaleReconstruct3D, "$Revision: 1.7 $");
vtkStandardNewMacro(vtkImageGreyscaleReconstruct3D);

struct coordAndOffset
{
  int x;
  int y;
  int z;
  unsigned long offset;
};


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
  // first check that we can actually use inDatas
  // We can't check inDatas directly; the array is only valid up to
  // NumberOfInputs.  Checking this->GetInput(idx) is the correct way.
  if (this->GetInput(1))
    {
    inDatas[1]->UpdateInformation();
    // copy all metadata (I hope) from the J marker image
    outData->CopyStructure(inDatas[1]);
    }
  else
    {
    vtkErrorMacro(<< "Marker input not set.");
    }
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
template <class T>
void fillNpPtrs(T* pPtr, int *Npi,
                int x, int y, int z,
                int xdim, int ydim, int zdim,
                T** NpPtrs)
{
  for (int i = 0; i < 13; i++)
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
}


template <class T>
void fillNmPtrs(T* pPtr, int *Nmi,
                int x, int y, int z,
                int xdim, int ydim, int zdim,
                T** NmPtrs)
{
  for (int i = 0; i < 13; i++)
    {
    // calculate modified pointer
    NmPtrs[i] = pPtr + Nmi[i];
    }

  // now zero things which are over the boundary
  if (x == 0)
    {
    NmPtrs[2] = NmPtrs[5] = NmPtrs[8] = NmPtrs[11] = (T*)NULL;
    }
  else if (x == xdim - 1)
    {
    NmPtrs[0] = NmPtrs[3] = NmPtrs[6] = NmPtrs[9] =
      NmPtrs[12] = (T*)NULL;
    }
  if (y == 0)
    {
    NmPtrs[6] = NmPtrs[7] = NmPtrs[8] = (T*)NULL;
    }
  else if (y == ydim - 1)
    {
    NmPtrs[0] = NmPtrs[1] = NmPtrs[2] =
      NmPtrs[9] = NmPtrs[10] = NmPtrs[11] = (T*)NULL;
    }
  if (z == zdim - 1)
    {
    NmPtrs[0] = NmPtrs[1] = NmPtrs[2] = \
      NmPtrs[3] = NmPtrs[4] = NmPtrs[5] =                \
      NmPtrs[6] = NmPtrs[7] = NmPtrs[8] = (T*)NULL;
    }
  // there are no z-1 cases
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
  int NpiO[13][3] = {{-1, -1, -1},
                     { 0, -1, -1},
                     {+1, -1, -1},
                     {-1,  0, -1},
                     { 0,  0, -1},
                     {+1,  0, -1},
                     {-1, +1, -1},
                     { 0, +1, -1},
                     {+1, +1, -1},
                     {-1, -1,  0},
                     { 0, -1,  0},
                     {+1, -1,  0},
                     {-1,  0,  0}};

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

  // Nminus neighbourhood pointers will be stored in this
  T *NmPtrs[13];
  // these are the Nmi offsets we need to extract the neighbourhood
  int Nmi[13];
  int NmiO[13][3] = {{+1, +1, +1},
                     { 0, +1, +1},
                     {-1, +1, +1},
                     {+1,  0, +1},
                     { 0,  0, +1},
                     {-1,  0, +1},
                     {+1, -1, +1},
                     { 0, -1, +1},
                     {-1, -1, +1},
                     {+1, +1,  0},
                     { 0, +1,  0},
                     {-1, +1,  0},
                     {+1,  0,  0}};
  

  Nmi[0] = +ydim*xdim + xdim + 1; // x+1, y+1, z+1
  Nmi[1] = +ydim*xdim + xdim;     // x,   y+1, z+1
  Nmi[2] = +ydim*xdim + xdim - 1; // x-1, y+1, z+1

  Nmi[3] = +ydim*xdim + 1;        // x+1, y,   z+1
  Nmi[4] = +ydim*xdim;            // x,   y,   z+1
  Nmi[5] = +ydim*xdim - 1;        // x-1, y,   z+1

  Nmi[6] = +ydim*xdim - xdim + 1; // x+1, y-1, z+1
  Nmi[7] = +ydim*xdim - xdim;     // x,   y-1, z+1
  Nmi[8] = +ydim*xdim - xdim - 1; // x-1, y-1, z+1

  Nmi[9] = +xdim + 1;             // x+1, y+1, z
  Nmi[10] = +xdim;                // x,   y+1, z
  Nmi[11] = +xdim - 1;            // x-1, y+1, z

  Nmi[12] = +1;                   // x+1, y,   z
  

  T *pPtr;
  T *Iptr;
  T minNp;
  T maxNp;
  T minNm;
  T maxNm;
  int x,y,z,i;
  // fifo of long offsets to the pixel in question
  vtkstd::queue<coordAndOffset> fifo;
  coordAndOffset tempCoordAndOffset;

  float progress = 0.0;
  
  if (Dual)
    {
    // 1. scan D1 in raster order
    pPtr = (T*)(outData->GetScalarPointer());
    Iptr = in1Ptr;
    for (z = 0; z < zdim; z++)
      {
      for (y = 0; y < ydim; y++)
        {
        for (x = 0; x < xdim; x++)
          {
          // J(p) max(min(J(q), q E N+(p) U p) ,I(p))
          // i.e. find the minimum in the N+ neighbourhood of P and P
          // replace J(p = x,y,z) with the max of that minimum and
          // I(p)

          // fill out structure with pointers to neighbours... if the
          // neighbour doesn't exist, the pointer will be NULL
          fillNpPtrs((T*)pPtr, Npi, x, y, z, xdim, ydim, zdim, (T**)NpPtrs);

          // a. determine minimum of Np neighbourhood
          minNp = *pPtr; // p itself is also a candidate
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
      // once every z increment, update the progress
      self->UpdateProgress((float)z / (float)(zdim - 1) * 0.33);
      }
    
    // 2. scan D1 in anti-raster order
    unsigned long pOffset = xdim * ydim * zdim - 1;    
    pPtr = (T*)(outData->GetScalarPointer()) + pOffset;
    Iptr = in1Ptr + pOffset;

    for (z = zdim - 1; z >= 0; z--)
      {
      for (y = ydim - 1; y >= 0; y--)
        {
        for (x = xdim - 1; x >= 0; x--)
          {
          // J(p) max(min(J(q), q E N-(p) U p) ,I(p))
          // i.e. find the minimum in the N- neighbourhood of P and P
          // replace J(p = x,y,z) with the max of that minimum and
          // I(p)

          // fill out NmPtrs for this position
          fillNmPtrs((T*)pPtr, Nmi, x, y, z, xdim, ydim, zdim, (T**)NmPtrs);

          // a. determine minimum of Np neighbourhood
          // now determine the minimum
          minNm = *pPtr; // p itself also counts
          for (i = 0; i < 13; i++)
            {
            if (NmPtrs[i] && *NmPtrs[i] < minNm)
              {
              minNm = *NmPtrs[i];
              }
            }
          *pPtr = minNm > *Iptr ? minNm : *Iptr;

          // extra step with anti-raster traversal: we have to start
          // filling the fifo...
          // search for J(q) in N-(p) such that J(q) > J(p) and
          // J(q) > I(q); if such a q exists, store *p* in the fifo
          for (i = 0; i < 13; i++)
            {
            if (NmPtrs[i] &&
                *NmPtrs[i] > *pPtr &&
                *NmPtrs[i] > *(Iptr + Nmi[i]))
              {
              tempCoordAndOffset.x = x;
              tempCoordAndOffset.y = y;
              tempCoordAndOffset.z = z;
              tempCoordAndOffset.offset = pOffset;
              fifo.push(tempCoordAndOffset);
              break;
              }
            }
          
          
          // decrement pointers
          pPtr--;
          Iptr--;
          // decrement offset
          pOffset--;
          }
        }
      // eventually two-thirds complete
      self->UpdateProgress(0.33 + (float)(zdim - 1 - z) / (float)(zdim - 1) * 0.33);
      }
    

    // 3. propagation step
    while (! fifo.empty())
      {
      // read the value from the fifo
      tempCoordAndOffset = fifo.front();
      // and take it off the queue
      fifo.pop();

      // set up our pointers
      pPtr = (T*)(outData->GetScalarPointer()) + tempCoordAndOffset.offset;
      Iptr = in1Ptr + tempCoordAndOffset.offset;
      x = tempCoordAndOffset.x;
      y = tempCoordAndOffset.y;
      z = tempCoordAndOffset.z;
      long currentOffset = tempCoordAndOffset.offset;

      // let's generate the complete neighbourhood, sanity checks and
      // all
      fillNpPtrs((T*)pPtr, Npi, x, y, z, xdim, ydim, zdim, (T**)NpPtrs);
      fillNmPtrs((T*)pPtr, Nmi, x, y, z, xdim, ydim, zdim, (T**)NmPtrs);

      for (i = 0; i < 13; i++)
        {
        // first for Np
        if (NpPtrs[i] &&
            *NpPtrs[i] > *pPtr &&
            *(Iptr + Npi[i]) !=  *NpPtrs[i])
          {
          *NpPtrs[i] = *pPtr > *(Iptr + Npi[i]) ? *pPtr : *(Iptr + Npi[i]);
          // current offset is that of "p", by adding Npi[i] we get
          // the offset of q
          tempCoordAndOffset.offset = currentOffset + Npi[i];
          tempCoordAndOffset.x = x + NpiO[i][0];
          tempCoordAndOffset.y = y + NpiO[i][1];
          tempCoordAndOffset.z = z + NpiO[i][2];
          fifo.push(tempCoordAndOffset);
          }
        // then for Nm
        if (NmPtrs[i] &&
            *NmPtrs[i] > *pPtr &&
            *(Iptr + Nmi[i]) !=  *NmPtrs[i])
          {
          *NmPtrs[i] = *pPtr > *(Iptr + Nmi[i]) ? *pPtr : *(Iptr + Nmi[i]);
          // current offset is that of "p", by adding Npi[i] we get
          // the offset of q
          tempCoordAndOffset.offset = currentOffset + Nmi[i];
          tempCoordAndOffset.x = x + NmiO[i][0];
          tempCoordAndOffset.y = y + NmiO[i][1];
          tempCoordAndOffset.z = z + NmiO[i][2];          
          fifo.push(tempCoordAndOffset);
          }
        }
      } // while (!fifo.empty() ...

    } // if (Dual) ...
  else
    {

    // ----------------------------------------------------------------------
    // Fast Hybrid Greyscale Reconstruction ---------------------------------
    // ----------------------------------------------------------------------

    // 1. scan D1 in raster order
    pPtr = (T*)(outData->GetScalarPointer());
    Iptr = in1Ptr;
    for (z = 0; z < zdim; z++)
      {
      for (y = 0; y < ydim; y++)
        {
        for (x = 0; x < xdim; x++)
          {
          // J(p) min(max(J(q), q E N+(p) U p), I(p))
          // i.e. find the maximum in the N+ neighbourhood of P and P
          // replace J(p = x,y,z) with the min of that maximum and
          // I(p)

          // fill out structure with pointers to neighbours... if the
          // neighbour doesn't exist, the pointer will be NULL
          fillNpPtrs((T*)pPtr, Npi, x, y, z, xdim, ydim, zdim, (T**)NpPtrs);

          // a. determine maximum of Np neighbourhood
          maxNp = *pPtr; // p itself is also a candidate
          for (i = 0; i < 13; i++)
            {
            if (NpPtrs[i] && *NpPtrs[i] > maxNp)
              {
              maxNp = *NpPtrs[i];
              }
            }
          *pPtr = maxNp < *Iptr ? maxNp : *Iptr;
          
          // increment pointer
          pPtr++;
          Iptr++;
          }
        }
      // once every z increment, update the progress
      self->UpdateProgress((float)z / (float)(zdim - 1) * 0.33);
      }
    
    // 2. scan D1 in anti-raster order
    unsigned long pOffset = xdim * ydim * zdim - 1;    
    pPtr = (T*)(outData->GetScalarPointer()) + pOffset;
    Iptr = in1Ptr + pOffset;

    for (z = zdim - 1; z >= 0; z--)
      {
      for (y = ydim - 1; y >= 0; y--)
        {
        for (x = xdim - 1; x >= 0; x--)
          {
          // J(p) <- min(max(J(q), q E N-(p) U p) ,I(p))
          // i.e. find the maximum in the N- neighbourhood of P and P
          // replace J(p = x,y,z) with the min of that maximum and
          // I(p)

          // fill out NmPtrs for this position
          fillNmPtrs((T*)pPtr, Nmi, x, y, z, xdim, ydim, zdim, (T**)NmPtrs);

          // a. determine minimum of Np neighbourhood
          // now determine the minimum
          maxNm = *pPtr; // p itself also counts
          for (i = 0; i < 13; i++)
            {
            if (NmPtrs[i] && *NmPtrs[i] > maxNm)
              {
              maxNm = *NmPtrs[i];
              }
            }
          *pPtr = maxNm < *Iptr ? maxNm : *Iptr;

          // extra step with anti-raster traversal: we have to start
          // filling the fifo...
          // search for J(q) in N-(p) such that J(q) > J(p) and
          // J(q) > I(q); if such a q exists, store *p* in the fifo
          for (i = 0; i < 13; i++)
            {
            if (NmPtrs[i] &&
                *NmPtrs[i] < *pPtr &&
                *NmPtrs[i] < *(Iptr + Nmi[i]))
              {
              tempCoordAndOffset.x = x;
              tempCoordAndOffset.y = y;
              tempCoordAndOffset.z = z;
              tempCoordAndOffset.offset = pOffset;
              fifo.push(tempCoordAndOffset);
              break;
              }
            }
          
          
          // decrement pointers
          pPtr--;
          Iptr--;
          // decrement offset
          pOffset--;
          }
        }
      // eventually two-thirds complete
      self->UpdateProgress(0.33 + (float)(zdim - 1 - z) / (float)(zdim - 1) * 0.33);
      }
    

    // 3. propagation step
    while (! fifo.empty())
      {
      // read the value from the fifo
      tempCoordAndOffset = fifo.front();
      // and take it off the queue
      fifo.pop();

      // set up our pointers
      pPtr = (T*)(outData->GetScalarPointer()) + tempCoordAndOffset.offset;
      Iptr = in1Ptr + tempCoordAndOffset.offset;
      x = tempCoordAndOffset.x;
      y = tempCoordAndOffset.y;
      z = tempCoordAndOffset.z;
      long currentOffset = tempCoordAndOffset.offset;

      // let's generate the complete neighbourhood, sanity checks and
      // all
      fillNpPtrs((T*)pPtr, Npi, x, y, z, xdim, ydim, zdim, (T**)NpPtrs);
      fillNmPtrs((T*)pPtr, Nmi, x, y, z, xdim, ydim, zdim, (T**)NmPtrs);

      for (i = 0; i < 13; i++)
        {
        // first for Np
        if (NpPtrs[i] &&
            *NpPtrs[i] < *pPtr &&
            *(Iptr + Npi[i]) !=  *NpPtrs[i])
          {
          *NpPtrs[i] = *pPtr < *(Iptr + Npi[i]) ? *pPtr : *(Iptr + Npi[i]);
          // current offset is that of "p", by adding Npi[i] we get
          // the offset of q
          tempCoordAndOffset.offset = currentOffset + Npi[i];
          tempCoordAndOffset.x = x + NpiO[i][0];
          tempCoordAndOffset.y = y + NpiO[i][1];
          tempCoordAndOffset.z = z + NpiO[i][2];
          fifo.push(tempCoordAndOffset);
          }
        // then for Nm
        if (NmPtrs[i] &&
            *NmPtrs[i] < *pPtr &&
            *(Iptr + Nmi[i]) !=  *NmPtrs[i])
          {
          *NmPtrs[i] = *pPtr < *(Iptr + Nmi[i]) ? *pPtr : *(Iptr + Nmi[i]);
          // current offset is that of "p", by adding Npi[i] we get
          // the offset of q
          tempCoordAndOffset.offset = currentOffset + Nmi[i];
          tempCoordAndOffset.x = x + NmiO[i][0];
          tempCoordAndOffset.y = y + NmiO[i][1];
          tempCoordAndOffset.z = z + NmiO[i][2];          
          fifo.push(tempCoordAndOffset);
          }
        }
      } // while (!fifo.empty() ...
    
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



