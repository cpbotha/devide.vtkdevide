#ifndef __vtkSelectConnectedComponents_h
#define __vtkSelectConnectedComponents_h

#include "vtkImageToImageFilter.h"
#include "vtkdevideImagingWin32Header.h"

class vtkImageConnector;
class vtkImageConnectorSeed;

// .NAME vtkImageSeedConnectivity - SeedConnectivity with user defined seeds.
// .SECTION Description
// vtkSelectedConnectedComponents marks pixels to user supplied seeds.
// The input must be unsigned long and the output will be unsigned
// char.  This filter is to be used when connected components (for
// example the result of a watershedding) are selected by setting seed
// points.  All selected connected components will be changed to
// OutputConnectedValue, all others to OutputUnconnectedValue.
// .AUTHOR Charl P. Botha
class VTK_DEVIDE_IMAGING_EXPORT vtkSelectConnectedComponents :
  public vtkImageToImageFilter
{
public:
  static vtkSelectConnectedComponents *New();
  vtkTypeRevisionMacro(vtkSelectConnectedComponents,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Methods for manipulating the seed pixels.
  void RemoveAllSeeds();
  void AddSeed(int num, int *index);
  void AddSeed(int i0, int i1, int i2);
  void AddSeed(int i0, int i1);

  // Description:
  // Set/Get what value is considered as connecting pixels.
  //vtkSetMacro(InputConnectValue, int);
  //vtkGetMacro(InputConnectValue, int);

  // Description:
  // Set/Get the value to set connected pixels to.
  vtkSetMacro(OutputConnectedValue, int);
  vtkGetMacro(OutputConnectedValue, int);

  // Description:
  // Set/Get the value to set unconnected pixels to.
  vtkSetMacro(OutputUnconnectedValue, int);
  vtkGetMacro(OutputUnconnectedValue, int);
  
  // Description:
  // Get the vtkImageCOnnector used by this filter.
  vtkGetObjectMacro(Connector,vtkImageConnector);

  // Description:
  // Set the number of axes to use in connectivity.
  vtkSetMacro(Dimensionality,int);
  vtkGetMacro(Dimensionality,int);
  
protected:
  vtkSelectConnectedComponents();
  ~vtkSelectConnectedComponents();

  //unsigned long InputConnectValue;
  unsigned char OutputConnectedValue;
  unsigned char OutputUnconnectedValue;
  vtkImageConnectorSeed *Seeds;
  vtkImageConnector *Connector;
  int Dimensionality;
  
  void ComputeInputUpdateExtents(vtkDataObject *out);

  void ExecuteInformation(vtkImageData *inData, vtkImageData *outData);
  void ExecuteData(vtkDataObject *out); 
private:
  vtkSelectConnectedComponents(const vtkSelectConnectedComponents&);  // Not implemented.
  void operator=(const vtkSelectConnectedComponents&);  // Not implemented.
};



#endif


  
