#include "vtkSelectConnectedComponents.h"

#include "vtkImageConnector.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"

vtkCxxRevisionMacro(vtkSelectConnectedComponents, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkSelectConnectedComponents);

//----------------------------------------------------------------------------
vtkSelectConnectedComponents::vtkSelectConnectedComponents()
{
//  this->InputConnectValue = 255;
  this->OutputConnectedValue = 255;
  this->OutputUnconnectedValue = 0;
  this->Seeds = NULL;
  this->Connector = vtkImageConnector::New();
  this->Dimensionality = 3;
}

//----------------------------------------------------------------------------
vtkSelectConnectedComponents::~vtkSelectConnectedComponents()
{
  this->Connector->Delete();
  this->RemoveAllSeeds();
}

//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::RemoveAllSeeds()
{
  vtkImageConnectorSeed *temp;
  while (this->Seeds)
    {
    temp = this->Seeds;
    this->Seeds = temp->Next;
    delete temp;
    }
}

//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::AddSeed(int num, int *index)
{
  int idx, newIndex[3];
  vtkImageConnectorSeed *seed;
  
  if (num > 3)
    {
    num = 3;
    } 
  for (idx = 0; idx < num; ++idx)
    {
    newIndex[idx] = index[idx];
    }
  for (idx = num; idx < 3; ++idx)
    {
    newIndex[idx] = 0;
    }
  seed = this->Connector->NewSeed(newIndex, NULL);
  seed->Next = this->Seeds;
  this->Seeds = seed;
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::AddSeed(int i0, int i1, int i2)
{
  int index[3];

  index[0] = i0;
  index[1] = i1;
  index[2] = i2;
  this->AddSeed(3, index);
}
//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::AddSeed(int i0, int i1)
{
  int index[2];

  index[0] = i0;
  index[1] = i1;
  this->AddSeed(2, index);
}

//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::ComputeInputUpdateExtents(vtkDataObject *out)
{
  vtkImageData *input = this->GetInput();

  out = out;
  if (input)
    {
    input->SetUpdateExtent(input->GetWholeExtent());
    }
}
//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::ExecuteInformation(
  vtkImageData *inData, vtkImageData *outData)
{
  // the default ExecuteInformation sets the same type on the output!
  outData->SetScalarType(VTK_UNSIGNED_CHAR);
}

//----------------------------------------------------------------------------
void vtkSelectConnectedComponents::ExecuteData(vtkDataObject *)
{
  vtkImageData *inData = this->GetInput();
  vtkImageData *outData = this->GetOutput();
  vtkImageConnectorSeed *seed;
  int idx0, idx1, idx2;
  vtkIdType inInc0, inInc1, inInc2;
  vtkIdType outInc0, outInc1, outInc2;
  int min0, max0, min1, max1, min2, max2;
  unsigned long *inPtr0, *inPtr1, *inPtr2;
  unsigned char *outPtr0, *outPtr1, *outPtr2;
  unsigned char temp1, temp2;
  int temp;

  outData->SetExtent(this->GetOutput()->GetWholeExtent());
  outData->AllocateScalars();
  
  if (inData->GetScalarType() != VTK_UNSIGNED_LONG ||
      outData->GetScalarType() != VTK_UNSIGNED_CHAR)
    {
    vtkErrorMacro("Execute: Input has to be unsigned long, output has to be "
                  "unsigned char.");
    return;
    }

  // Pick an intermediate value (In some cases, we could eliminate the last threshold.)
  temp1 = 1;
  while (temp1 == this->OutputUnconnectedValue ||
         temp1 == this->OutputConnectedValue)
    {
    ++temp1;
    }
  temp2 = temp1 + 1;
  while (temp2 == this->OutputUnconnectedValue ||
         temp2 == this->OutputConnectedValue)
    {
    ++temp2;
    }

  // -----
  // See how many seeds there are and extract the scalars at their
  // positions...
  int numSeeds = 0;
  vtkImageConnectorSeed *tempSeed;
  
  tempSeed = this->Seeds;
  while (tempSeed)
    {
    ++numSeeds;
    tempSeed = tempSeed->Next;
    }

  unsigned long *inputConnectedValues = new unsigned long [numSeeds];

  this->GetInput()->GetExtent(min0, max0, min1, max1, min2, max2);
  tempSeed = this->Seeds;
  for (int seedIdx = 0; seedIdx < numSeeds; seedIdx++)
    {

    // make sure z value of seed is acceptable
    if (tempSeed->Index[2] < min2)
      {
      tempSeed->Index[2] = min2;
      }
    if (tempSeed->Index[2] > max2)
      {
      tempSeed->Index[2] = max2;
      }

    // get scalar value at that point
    inputConnectedValues[seedIdx] = 
      *((unsigned long *)(inData->GetScalarPointer(tempSeed->Index)));

    // and go to the next seed for the next iteration
    tempSeed = tempSeed->Next;
    }

  // now we've built up a list with scalar values at the seed positions
  // -----
  
  
  
  //-------
  // threshold to eliminate unknown values ( only intermediate and 0)
  inData->GetIncrements(inInc0, inInc1, inInc2);
  this->GetOutput()->GetExtent(min0, max0, min1, max1, min2, max2);
  outData->GetIncrements(outInc0, outInc1, outInc2);
  inPtr2 = (unsigned long *)(inData->GetScalarPointer(min0,min1,min2));
  outPtr2 = (unsigned char *)(outData->GetScalarPointer(min0,min1,min2));

  int tempSeedIdx;
  int valueIsSeedValue;
  
  for (idx2 = min2; idx2 <= max2; ++idx2)
    {
    inPtr1 = inPtr2;
    outPtr1 = outPtr2;
    for (idx1 = min1; idx1 <= max1; ++idx1)
      {
      inPtr0 = inPtr1;
      outPtr0 = outPtr1;
      for (idx0 = min0; idx0 <= max0; ++idx0)
        {
        
        // check if *inPtr0 is equal to ONE of the seedvalues
        // instead of the old inputConnectedValue:
        // "if (*inPtr0 == this->InputConnectValue)"
        
        valueIsSeedValue = 0;
        for (tempSeedIdx = 0; tempSeedIdx < numSeeds && !valueIsSeedValue;
             tempSeedIdx++)
          {
          if (inputConnectedValues[tempSeedIdx] == *inPtr0)
            {
            valueIsSeedValue = 1;
            }
          }
        
        if (valueIsSeedValue)
          {
          *outPtr0 = temp1;
          }
        else
          {
          *outPtr0 = 0;
          }
        
        inPtr0 += inInc0;
        outPtr0 += outInc0;
        }
      inPtr1 += inInc1;
      outPtr1 += outInc1;
      }
    inPtr2 += inInc2;
    outPtr2 += outInc2;
    }

  // get rid of this heap allocation before we possibly abort
  delete[] inputConnectedValues;
  
  this->UpdateProgress(0.2);
  if (this->AbortExecute)
    {
    return;
    }
  
  //-------
  // find actual seeds in this image. (only scan along the first axis for now)
  this->Connector->RemoveAllSeeds();
  seed = this->Seeds;
  while (seed)
    {
    temp = seed->Index[0];
    // make sure z value of seed is acceptable
    if (seed->Index[2] < min2)
      {
      seed->Index[2] = min2;
      }
    if (seed->Index[2] > max2)
      {
      seed->Index[2] = max2;
      }
    outPtr0 = (unsigned char *)(outData->GetScalarPointer(seed->Index));
    for (idx0 = temp; idx0 <= max0; ++idx0)
      {
      if (*outPtr0 == temp1)
        { // we found our seed
        seed->Index[0] = idx0;
        this->Connector->AddSeed(this->Connector->NewSeed(seed->Index, outPtr0));
        seed->Index[0] = temp;
        break;
        }
      outPtr0 += outInc0;
      }
    seed = seed->Next;
    }

  this->UpdateProgress(0.5);
  if (this->AbortExecute)
    {
    return;
    }

  //-------
  // connect
  this->Connector->SetUnconnectedValue(temp1);
  this->Connector->SetConnectedValue(temp2);
  this->Connector->MarkData(outData, this->Dimensionality, 
                            this->GetOutput()->GetExtent());

  this->UpdateProgress(0.9);
  if (this->AbortExecute)
    {
    return;
    }

  //-------
  // Threshold to convert intermediate values into OutputUnconnectedValues
  outPtr2 = (unsigned char *)(outData->GetScalarPointer(min0,min1,min2));
  for (idx2 = min2; idx2 <= max2; ++idx2)
    {
    outPtr1 = outPtr2;
    for (idx1 = min1; idx1 <= max1; ++idx1)
      {
      outPtr0 = outPtr1;
      for (idx0 = min0; idx0 <= max0; ++idx0)
        {
        if (*outPtr0 == temp2)
          {
          *outPtr0 = this->OutputConnectedValue;
          }
        else
          {
          *outPtr0 = this->OutputUnconnectedValue;
          }
        outPtr0 += outInc0;
        }
      outPtr1 += outInc1;
      }
     outPtr2 += outInc2;
    }
}


void vtkSelectConnectedComponents::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->Connector )
    {
    os << indent << "Connector: " << this->Connector << "\n";
    }
  else
    {
    os << indent << "Connector: (none)\n";
    }

  os << indent << "Dimensionality: " << this->Dimensionality << "\n";
//  os << indent << "InputConnectValue: " << this->InputConnectValue << "\n";
  os << indent << "OutputConnectedValue: " << this->OutputConnectedValue << "\n";
  os << indent << "OutputUnconnectedValue: " << this->OutputUnconnectedValue << "\n";
}

