#ifndef __vtkHistogramLookupTable_h
#define __vtkHistogramLookupTable_h

#include "vtkImageTwoInputFilter.h"
#include "vtkdevideImagingWin32Header.h"


class VTK_DEVIDE_IMAGING_EXPORT vtkHistogramLookupTable :
public vtkImageTwoInputFilter
{
public:
  static vtkHistogramLookupTable *New();
  vtkTypeRevisionMacro(vtkHistogramLookupTable,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  
protected:
  vtkHistogramLookupTable();
  ~vtkHistogramLookupTable() {};
  
  void ExecuteInformation(vtkImageData **inDatas, vtkImageData *outData);
  virtual void ComputeInputUpdateExtent(int inExt[6], int outExt[6],
                                        int whichInput);
  // this is called by the SuperClass - we just chain to the other EI
  void ExecuteInformation(){this->vtkImageTwoInputFilter::ExecuteInformation();};
  
private:
  void ExecuteData(vtkDataObject *out);
  
  vtkHistogramLookupTable(const vtkHistogramLookupTable&);  // Not implemented.
  void operator=(const vtkHistogramLookupTable&);  // Not implemented.
};

#endif



