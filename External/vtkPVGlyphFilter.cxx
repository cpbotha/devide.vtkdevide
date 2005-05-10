/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkPVGlyphFilter.cxx,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPVGlyphFilter.h"

#include "vtkGarbageCollector.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMaskPoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkCxxRevisionMacro(vtkPVGlyphFilter, "$Revision: 1.10 $");
vtkStandardNewMacro(vtkPVGlyphFilter);

//-----------------------------------------------------------------------------
vtkPVGlyphFilter::vtkPVGlyphFilter()
{
  this->SetColorModeToColorByScalar();
  this->SetScaleModeToScaleByVector();
  this->MaskPoints = vtkMaskPoints::New();
  this->MaximumNumberOfPoints = 5000;
  this->UseMaskPoints = 1;
}

//-----------------------------------------------------------------------------
vtkPVGlyphFilter::~vtkPVGlyphFilter()
{
  if(this->MaskPoints)
    {
    this->MaskPoints->Delete();
    }
}

//-----------------------------------------------------------------------------
void vtkPVGlyphFilter::SetRandomMode(int mode)
{
  this->MaskPoints->SetRandomMode(mode);
}

//-----------------------------------------------------------------------------
int vtkPVGlyphFilter::GetRandomMode()
{
  return this->MaskPoints->GetRandomMode();
}

//-----------------------------------------------------------------------------
int vtkPVGlyphFilter::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  if (this->UseMaskPoints)
    {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet* input = vtkDataSet::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet* inputCopy = input->NewInstance();
    inputCopy->ShallowCopy(input);
    this->MaskPoints->SetInput(inputCopy);
    inputCopy->Delete();

    vtkIdType maxNumPts = this->MaximumNumberOfPoints;
    vtkIdType numPts = inputCopy->GetNumberOfPoints();

    // Although this is not perfectly process invariant, it is better
    // than we had before (divide by number of processes).
    vtkIdType totalNumPts = numPts;

    maxNumPts = (maxNumPts < 1) ? 1 : maxNumPts;
    this->MaskPoints->SetMaximumNumberOfPoints(maxNumPts);
    this->MaskPoints->SetOnRatio(numPts / maxNumPts);
    // I do not like connecting internal filters to the actual input, but
    // This is the smallest change possible to fix the problem.
    // This update caused input to be executed with number of pieces of 1.
    vtkInformation *maskPointsInfo =
      this->MaskPoints->GetOutputPortInformation(0);
    maskPointsInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
    maskPointsInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
    maskPointsInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
    this->MaskPoints->Update();
    }

  vtkInformationVector* inputVs[2];

  vtkInformationVector* inputV = inputVector[0];
  inputVs[0] = vtkInformationVector::New();
  inputVs[0]->SetNumberOfInformationObjects(1);
  vtkInformation* inInfo = vtkInformation::New();
  inInfo->Copy(inputV->GetInformationObject(0));
  inInfo->Set(vtkDataObject::DATA_OBJECT(), this->MaskPoints->GetOutput());
  inputVs[0]->SetInformationObject(0, inInfo);
  inInfo->Delete();
  inputVs[1] = inputVector[1];
  
  // added by cpbotha:
  // if the vector selection is at the default of NULL and there are no 
  // vectors in the input dataset, see if there are any scalars.  If
  // so, select these as the input vectors (they could have enough
  // elements to be handled as vectors)

  // funky new-style way of finding input
//  vtkDataSet* input = vtkDataSet::SafeDownCast(
//    inInfo->Get(vtkDataObject::DATA_OBJECT()));
//
//  if (input)
//    {
//    vtkPointData *pd = input->GetPointData();
//    if (!pd->GetVectors(this->InputVectorsSelection) &&
//        !this->InputVectorsSelection)
//      {
      // NO vectors, so let's try selecting the scalars (if there
      // are any!)
//      vtkDataArray *inScalars = pd->GetScalars(
//        this->InputScalarsSelection);
//
//      if (inScalars)
//        {
//        this->SelectInputVectors(inScalars->GetName());
//        }
//      }
//    }
  // end addition by cpbotha
  
   
  
  int retVal = this->Superclass::RequestData(request, inputVs, outputVector);
  inputVs[0]->Delete();

  return retVal;
}

//-----------------------------------------------------------------------------
void vtkPVGlyphFilter::ReportReferences(vtkGarbageCollector* collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->MaskPoints, "MaskPoints");
}

//-----------------------------------------------------------------------------
void vtkPVGlyphFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  
  os << indent << "MaximumNumberOfPoints: " << this->GetMaximumNumberOfPoints()
     << endl;

  os << indent << "UseMaskPoints: " << (this->UseMaskPoints?"on":"off") << endl;

}
