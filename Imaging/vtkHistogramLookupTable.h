#ifndef __vtkHistogramLookupTable_h
#define __vtkHistogramLookupTable_h

#include "vtkImageTwoInputFilter.h"
#include "vtkdevideImagingWin32Header.h"

// FIXME:
// add this to cvs
// we're going to take two inputs: a multi-component image and the
// histogram... the output will be a single-component volume


class VTK_DEVIDE_IMAGING_EXPORT vtkHistogramLookupTable :
public vtkImageTwoInputFilter
{
public:
  static vtkHistogramLookupTable *New();
  vtkTypeRevisionMacro(vtkHistogramLookupTable,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  vtkSetMacro(Input1Bins,int);
  vtkGetMacro(Input1Bins,int);

  vtkSetMacro(Input2Bins,int);
  vtkGetMacro(Input2Bins,int);

  vtkSetMacro(MaxSamplesPerBin,long);
  vtkGetMacro(MaxSamplesPerBin,long);

  vtkSetMacro(ClipSamples, int);
  vtkGetMacro(ClipSamples, int);
  vtkBooleanMacro(ClipSamples, int);
  
protected:
  vtkHistogramLookupTable();
  ~vtkHistogramLookupTable() {};

  int Input1Bins;
  int Input2Bins;

  long MaxSamplesPerBin;

  int ClipSamples;
  
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



