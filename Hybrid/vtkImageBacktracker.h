
#ifndef __vtkImageBacktracker_h
#define __vtkImageBacktracker_h

#include "vtkImageToPolyDataFilter.h"
#include "vtkdevideHybridWin32Header.h"

class vtkImageConnectorSeed;

//BTX

//ETX

class VTK_DEVIDE_HYBRID_EXPORT vtkImageBacktracker :
public vtkImageToPolyDataFilter
{
public:
  static vtkImageBacktracker *New();
  vtkTypeRevisionMacro(vtkImageBacktracker,vtkImageToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  void RemoveAllSeeds();
  void AddSeed(int num, int *index);
  void AddSeed(int i0, int i1, int i2);
  void AddSeed(int i0, int i1);

  //  bool in_range( int x, int y, int z ) { return x>

  vtkImageConnectorSeed *Seeds;

protected:
  vtkImageBacktracker();
  ~vtkImageBacktracker();


  void ExecuteInformation(vtkImageData *inData, vtkPolyData *outData);

  virtual void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);
  // this is called by the SuperClass - we just chain to the other EI
  void ExecuteInformation(){this->vtkImageToPolyDataFilter::ExecuteInformation();};
  
private:
  void ExecuteData(vtkDataObject *out);
  vtkImageConnectorSeed *NewSeed(int index[3], void *ptr);
  
  vtkImageBacktracker(const vtkImageBacktracker&);  // Not implemented.
  void operator=(const vtkImageBacktracker&);  // Not implemented.
};

#endif



