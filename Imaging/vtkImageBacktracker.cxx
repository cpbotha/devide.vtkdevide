#include "vtkImageBacktracker.h"
#include "vtkImageConnector.h"
#include "vtkImageData.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"

#include <vtkstd/list>
#include <vtkstd/map>
#include <vtkstd/vector>
#include <vtkstd/queue>


vtkCxxRevisionMacro(vtkImageBacktracker, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImageBacktracker);

//----------------------------------------------------------------------------
vtkImageBacktracker::vtkImageBacktracker()
{
    this->Seeds = NULL;
}

// ----
vtkImageBacktracker::~vtkImageBacktracker()
{
  this->RemoveAllSeeds();
}

//----------------------------------------------------------------------------
void vtkImageBacktracker::RemoveAllSeeds()
{
  vtkImageConnectorSeed *temp;
  while (this->Seeds)
    {
    temp = this->Seeds;
    this->Seeds = temp->Next;
    delete temp;
    }
}

// -------------------------------------------------
vtkImageConnectorSeed *vtkImageBacktracker::NewSeed(int index[3], void *ptr)
{
  vtkImageConnectorSeed *seed = vtkImageConnectorSeed::New();
  int idx;

  for (idx = 0; idx < 3; ++idx)
    {
    seed->Index[idx] = index[idx];
    }
  seed->Pointer = ptr;
  seed->Next = NULL;

  return seed;
}


//----------------------------------------------------------------------------
void vtkImageBacktracker::AddSeed(int num, int *index)
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
  seed = NewSeed(newIndex, NULL);

  seed->Next = this->Seeds;
  this->Seeds = seed;
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkImageBacktracker::AddSeed(int i0, int i1, int i2)
{
  int index[3];

  index[0] = i0;
  index[1] = i1;
  index[2] = i2;
  this->AddSeed(3, index);
}
//----------------------------------------------------------------------------
void vtkImageBacktracker::AddSeed(int i0, int i1)
{
  int index[2];

  index[0] = i0;
  index[1] = i1;
  this->AddSeed(2, index);
}


//----------------------------------------------------------------------------
// Grow the output image 
void vtkImageBacktracker::ExecuteInformation(
                    vtkImageData *inData, vtkPolyData *outData)
{
  //    outData->SetNumberOfScalarComponents(1);
  //    outData->SetScalarType(VTK_DOUBLE);

    //    inDatas[1]->GetScalarRange(in2range);

    double origin[3];
    inData->GetOrigin( origin );
    //    outData->SetOrigin( origin );
    double spacing[3];
    inData->GetSpacing( spacing );
    //    outData->SetSpacing( spacing );
    
    int extent[6];
    inData->GetExtent( extent );
    //    outData->SetExtent( extent );
    //    outData->SetWholeExtent( extent );

}

//----------------------------------------------------------------------------
void vtkImageBacktracker::ComputeInputUpdateExtent(int inExt[6], 
                                                   int outExt[6])
{
    // get the whole image for input 2
    memcpy(inExt,this->GetInput()->GetWholeExtent(),6*sizeof(int));
}


// direction conversion routines
unsigned char pxpypzToDirection( int px, int py, int pz )
{
  return (px+1)+3*(py+1)+9*(pz+1);
}

void directionTopxpypz( unsigned char c, int& px, int& py, int& pz )
{
  px = (c % 3) - 1;
  py = ((c/3)%3) -1;
  pz = ((c/9)%3) -1;
}


//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
// Handles the two input operations
template <class T>
void vtkImageBacktrackerExecute(vtkImageBacktracker *self,
                                vtkImageData *inData, T *inPtr,
                                vtkPolyData *outData, double *outPtr)
{
  int* dims = inData->GetDimensions();
		

  vtkPoints *points = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();

  double origin[3];
  inData->GetOrigin( origin );
  double spacing[3];
  inData->GetSpacing( spacing );

  int i;
  for (  vtkImageConnectorSeed * s = self->Seeds;s!=NULL; s=s->Next) 
    {

      int ix,iy,iz;
      ix = s->Index[0];
      iy = s->Index[1];
      iz = s->Index[2];

      
      int s1=-1,s2=-1;

      vtkIdType fiber[5000];
      int length=0;

      while(length<5000)
        {
          float x[3];
          x[0]=ix*spacing[0]+origin[0];
          x[1]=iy*spacing[1]+origin[1];
          x[2]=iz*spacing[2]+origin[2];
          //          s1 = s2;
          s2 = points->InsertNextPoint(x); // starting point
          fiber[length] = s2;
          length++;


          int dx,dy,dz;
          int inpos[3] = {ix,iy,iz};
          vtkIdType source_cell_id = inData->ComputePointId( inpos );
          //          printf("cell nr %d == %d\n", source_cell_id, ix+iy*dims[1]+iz*dims[1]*dims[2] );
          directionTopxpypz( inPtr[source_cell_id], dx,dy,dz );
          if (dx==0 && dy==0 && dz==0)
            break;
          //          printf("ix,iy,iz: %d %d %d   dx,dy,dz: %d %d %d\n", ix,iy,iz, dx,dy,dz );
          ix-=dx; iy-=dy; iz-=dz;
        }

      if (length>1)
            lines->InsertNextCell(length,fiber);

    }

  outData->SetPoints(points);
  outData->SetLines( lines );

  points->Delete();
  lines->Delete();

  // clear output data
  //memset(outPtr, 0, input1bins * input2bins * sizeof(double));

  // first we iterate over the seed pixels, and add them to the Known list

//   for( vtkImageConnectorSeed * s = self->Seeds; s; s=s->Next )
//     {
//       known.push_back( vtkFastMarchingPoint( s->Index[0], s->Index[1], s->Index[2], 0 ));
//     }

  //    long MaxSamplesPerBin = self->GetMaxSamplesPerBin();
  //   for (unsigned long i = 0; i < noe; i++)
  //     {
  //       if (i % noeProgressStep == 0)
  //         {
  // 	  self->UpdateProgress(progress);
  // 	  progress += progressStep;
  // 	  // progress will end up one step higher than 1.0, but that's
  // 	  // okay, because we won't be using it then.
  //         }
  
  //       usleep(1000);
  
  //     }
  
  self->UpdateProgress(1.0);
}


void vtkImageBacktracker::ExecuteData(vtkDataObject *out)
{
  // Make sure the Input has been set.
  if ( this->GetInput() == NULL)    
    {
    vtkErrorMacro(<< "ExecuteData: No input set.");
    return;
    }
  
  // Too many filters have floating point exceptions to execute
  // with empty input/ no request.
  if (this->UpdateExtentIsEmpty(out))
    {
    return;
    }

  // get metadata across
  this->ExecuteInformation();
  // make sure output has been allocated
  vtkPolyData *output = vtkPolyData::SafeDownCast(out);

  //  output->AllocateScalars();
  
  // now let's start with the actual work
  void *in1Ptr = this->GetInput()->GetScalarPointer();
  double *outPtr = NULL;//  i have no clue what this was for .. (double*)output->GetScalarPointer();

  this->GetInput()->Update();

  switch (this->GetInput()->GetScalarType())
  {
      vtkTemplateMacro(
          vtkImageBacktrackerExecute(
              this,
              this->GetInput(), static_cast<VTK_TT*>(in1Ptr),
              output, outPtr
              )
          );

    default:                              
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }


  
}


void vtkImageBacktracker::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

//   os << indent << "Input1Bins: " << this->Input1Bins << "\n";
//   os << indent << "Input2Bins: " << this->Input2Bins << "\n";  
  
}




