// .NAME vtkGdcmSerieHelper - VTK class to encapsulate gdcm::seriehelper
// .SECTION Description
// GDCM supplies vtkGdcmReader and vtkGdcmSerieWriter, but no VTK
// encapsulation of the seriehelper, so I slapped this together.
// copyright 2006 Charl P. Botha <http://cpbotha.net/>

#ifndef __vtkGdcmSerieHelper_h
#define __vtkGdcmSerieHelper_h

#include "gdcm.h"
//#include "vtkGdcmReader.h"
#include "vtkObject.h"

class VTK_EXPORT vtkGdcmSerieHelper : public vtkObject
{
public:
  // Methods from vtkObject
  vtkTypeRevisionMacro(vtkGdcmSerieHelper,vtkObject);
  // Description:
  // Print ObjectFactor to stream.
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  static vtkGdcmSerieHelper* New();

  void SetDirectory(char *);
//BTX
  gdcm::FileList *GetFileList(char *);
//ETX
  int GetNumberOfUIDs(void);
  char *GetUID(int);
//  void SetFileListOnReader(vtkGdcmReader *reader,
//                           char *serieUID);

protected:
  vtkGdcmSerieHelper();
  virtual ~vtkGdcmSerieHelper();

  gdcm::SerieHelper *helper;

  char temp_uid[8192];

private:
  vtkGdcmSerieHelper(const vtkGdcmSerieHelper&);  // Not implemented.
  void operator=(const vtkGdcmSerieHelper&);  // Not implemented.  
};


#endif __vtkGdcmSerieHelper_h
