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
#include "vtkMaskPoints.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"

vtkCxxRevisionMacro(vtkPVGlyphFilter, "$Revision: 1.7 $");
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
void vtkPVGlyphFilter::SetInput(vtkDataSet *input)
{
  this->MaskPoints->SetInput(input);
  this->Superclass::SetInput(this->MaskPoints->GetOutput());
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
void vtkPVGlyphFilter::Execute()
{
  if (this->UseMaskPoints)
    {
    vtkPolyData* output = this->GetOutput();
    this->Superclass::SetInput(this->MaskPoints->GetOutput());
    vtkIdType maxNumPts = this->MaximumNumberOfPoints;
    vtkIdType numPts = this->MaskPoints->GetInput()->GetNumberOfPoints();
    // Although this is not perfectly process invariant, it is better
    // than we had before (divide by number of processes).
    vtkIdType totalNumPts = numPts;

    maxNumPts = (maxNumPts < 1) ? 1 : maxNumPts;
    this->MaskPoints->SetMaximumNumberOfPoints(maxNumPts);
    this->MaskPoints->SetOnRatio(numPts / maxNumPts);
    // I do not like connecting internal filters to the actual input, but
    // This is the smallest change possible to fix the problem.
    // This update caused input to be executed with number of piecces of 1.
    this->MaskPoints->GetOutput()->SetUpdateNumberOfPieces(output->GetUpdateNumberOfPieces());
    this->MaskPoints->GetOutput()->SetUpdatePiece(output->GetUpdatePiece());
    this->MaskPoints->GetOutput()->SetUpdateGhostLevel(output->GetUpdateGhostLevel());
    this->MaskPoints->Update();
    }
  else
    {
    this->Superclass::SetInput(this->MaskPoints->GetInput());
    }

  // added by cpbotha:
  // if the vector selection is at the default of NULL and there are no 
  // vectors in the input dataset, see if there are any scalars.  If
  // so, select these as the input vectors (they could have enough
  // elements to be handled as vectors)
  vtkDataSet *input = this->GetInput();
  if (input)
    {
    vtkPointData *pd = input->GetPointData();
    if (!pd->GetVectors(this->InputVectorsSelection) && 
        !this->InputVectorsSelection)
      {
      // NO vectors, so let's try selecting the scalars (if there
      // are any!)
      vtkDataArray *inScalars = pd->GetScalars(
        this->InputScalarsSelection);
      
      if (inScalars)
        {
        this->SelectInputVectors(inScalars->GetName());
        }
      }
    }
  // end addition by cpbotha
   
  this->Superclass::Execute();
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
  os << indent << "InputScalarsSelection: " 
     << (this->InputScalarsSelection ? this->InputScalarsSelection : "(none)")
     << endl;

  os << indent << "InputVectorsSelection: " 
     << (this->InputVectorsSelection ? this->InputVectorsSelection : "(none)")
     << endl;

  os << indent << "InputNormalsSelection: " 
     << (this->InputNormalsSelection ? this->InputNormalsSelection : "(none)")
     << endl;
  
  os << indent << "MaximumNumberOfPoints: " << this->GetMaximumNumberOfPoints()
     << endl;

  os << indent << "UseMaskPoints: " << (this->UseMaskPoints?"on":"off") << endl;
}
