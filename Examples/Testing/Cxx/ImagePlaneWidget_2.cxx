#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkDataSet.h"
#include "vtkCommand.h"
#include "vtkImageActor.h"
#include "vtkImageData.h"
#include "vtkImageMapToColors.h"
#include "vtkImagePlaneWidget2.h"
#include "vtkImageReader.h"
#include "vtkLookupTable.h"
#include "vtkScalarsToColors.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume16Reader.h"
#include "vtkTexture.h"



int main( void )
{
/*
  vtkImageReader*  v16 = vtkImageReader::New();
    v16->SetFileName("c:\\Builder\\vtkMyData\\sagg_full4.pci");
    v16->SetDataScalarTypeToUnsignedShort();
    v16->SetDataByteOrderToLittleEndian();
    v16->SetFileDimensionality(3);
    v16->FileLowerLeftOff();
    v16->SetDataExtent(0,511,0,511,0,59);
    v16->SetDataSpacing (0.3125, 0.3125, 2.0);
    v16->Update();
*/
  vtkVolume16Reader* v16 =  vtkVolume16Reader::New();
    v16->SetDataDimensions( 64, 64);
    v16->GetOutput()->SetOrigin( 0.0, 0.0, 0.0 );
    v16->SetDataByteOrderToLittleEndian();
    v16->SetFilePrefix( "/home/cpbotha/build/VTKData/Data/headsq/quarter" );
    v16->SetImageRange( 1, 93);
    v16->SetDataSpacing( 0.32, 0.32, 0.15);
    v16->Update();

  vtkOutlineFilter* outline = vtkOutlineFilter::New();
    outline->SetInput(v16->GetOutput());

  vtkPolyDataMapper* outlineMapper = vtkPolyDataMapper::New();
    outlineMapper->SetInput(outline->GetOutput());

  vtkActor* outlineActor =  vtkActor::New();
    outlineActor->SetMapper(outlineMapper);

  vtkRenderer* ren1 = vtkRenderer::New();
  vtkRenderWindow* renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);

  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

  vtkCellPicker* picker = vtkCellPicker::New();
    picker->SetTolerance(0.005);
    picker->PickFromListOn();

  vtkImagePlaneWidget2* planeWidgetX = vtkImagePlaneWidget2::New();
    planeWidgetX->DisplayTextOn();
    planeWidgetX->TextureInterpolateOff();
    planeWidgetX->SetInput(v16->GetOutput());
    planeWidgetX->SetResliceInterpolateToNearestNeighbour();
    planeWidgetX->SetInteractor(iren);
    planeWidgetX->SetPlaneOrientationToXAxes();
    planeWidgetX->SetSliceIndex(32);
    planeWidgetX->SetPicker(picker);
    planeWidgetX->GetPlaneProperty()->SetColor(1,0,0);
    planeWidgetX->On();

  vtkImagePlaneWidget2* planeWidgetY = vtkImagePlaneWidget2::New();
    planeWidgetY->DisplayTextOn();
    planeWidgetY->TextureInterpolateOff();
    planeWidgetY->SetInput(v16->GetOutput());
    planeWidgetY->SetResliceInterpolateToNearestNeighbour();
    planeWidgetY->SetInteractor(iren);
    planeWidgetY->SetPlaneOrientationToYAxes();
    planeWidgetY->SetSliceIndex(32);
    planeWidgetY->SetPicker(picker);
    planeWidgetY->GetPlaneProperty()->SetColor(1,1,0);
    planeWidgetY->On();

  vtkImagePlaneWidget2* planeWidgetZ = vtkImagePlaneWidget2::New();
    planeWidgetZ->DisplayTextOn();
    planeWidgetZ->TextureInterpolateOff();
    planeWidgetZ->SetInput(v16->GetOutput());
    planeWidgetZ->SetResliceInterpolateToNearestNeighbour();
    planeWidgetZ->SetInteractor(iren);
    planeWidgetZ->SetPlaneOrientationToZAxes();
    planeWidgetZ->SetSliceIndex(46);
    planeWidgetZ->SetPicker(picker);
    planeWidgetZ->GetPlaneProperty()->SetColor(0,0,1);
    planeWidgetZ->On();


  ren1->AddActor( outlineActor);
  ren1->SetBackground( 0.1, 0.1, 0.2);

  renWin->SetSize(700,700);

  ren1->GetActiveCamera()->Elevation(110);
  ren1->GetActiveCamera()->SetViewUp(0, 0, -1);
  ren1->GetActiveCamera()->Azimuth(45);
  ren1->ResetCameraClippingRange();


  renWin->Render();
  // force a repostion of all actors in all renderers
  iren->SetEventPosition(350,350);
  iren->SetKeyCode('r');
  iren->InvokeEvent(vtkCommand::CharEvent,NULL);
  renWin->Render();

  iren->Start();

  planeWidgetX->Delete();
  planeWidgetY->Delete();
  planeWidgetZ->Delete();
  picker->Delete();
  outlineActor->Delete();
  outlineMapper->Delete();
  outline->Delete();
  iren->Delete();
  renWin->Delete();
  ren1->Delete();
  v16->Delete();

  return 0;
}


