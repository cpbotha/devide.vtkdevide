// vtkHDFVolumeReader.h -- Charl P. Botha <c.p.botha@its.tudelft.nl>
// $Id: vtkHDFVolumeReader.h,v 1.1 2003/01/08 14:07:29 cpbotha Exp $
// vtk class for reading scalar volume data

#ifndef vtkHDFVolumeReader_h
#define vtkHDFVolumeReader_h

#include <mfhdf.h>
#include <string>
#include <stdexcept>
#include <vtkStructuredPointsSource.h>

#include "vtkdscasIOWin32Header.h"

/**
 * Simple vtk class to read SDS HDF data.  It only reads the first SDS in
 * the HDF file.  At the moment, it assumes that the data is 3dimensional
 * and that it has been stored as unsigned shorts.
 * @author Charl P. Botha <c.p.botha@its.tudelft.nl>
 */
class VTK_DSCAS_IO_EXPORT vtkHDFVolumeReader : public vtkStructuredPointsSource
{
protected:
   float DataSpacing[3];
   float DataOrigin[3];
   int data_dimensions[3];
   char* FileName;
   int32 sd_id;
   int32 sds_id;
   
   vtkHDFVolumeReader();
   ~vtkHDFVolumeReader();
   /**
    * Opens hdf file and accesses first dataset.
    * 
    * The method calls close_hdf() before it gets to work, so it's idempotent.
    * It will open the hdf file and store file and dataset handles in sd_id
    * and sds_id.
    */
   int open_hdf(void);
   /// Idempotent function that makes sure any open HDF file is closed.
   void close_hdf(void);
   /**
    * Method required by the VTK API that reads and sets all metadata on
    * the output.
    */
   void ExecuteInformation();
   /**
    * Method required by the VTK API that reads the actual data.
    * 
    * By the time this method is called, ExecuteInformation will have been
    * called.  Please see the body of this method for details on how to work
    * with the passed parameter out.
    */
   void ExecuteData(vtkDataObject* out);
   
public:
   vtkTypeMacro(vtkHDFVolumeReader, vtkStructuredPointsSource);
   static vtkHDFVolumeReader *New();
   vtkSetStringMacro(FileName);
   vtkGetStringMacro(FileName);
   vtkSetVector3Macro(DataSpacing,float);
   vtkGetVectorMacro(DataSpacing,float,3);
   vtkSetVector3Macro(DataOrigin,float);
   vtkGetVectorMacro(DataOrigin,float,3);
};
#endif

