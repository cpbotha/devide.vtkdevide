// vtkDICOMVolumeReader.cxx copyright (c) 2001,2002 Charl P. Botha <cpbotha@ieee.org>
// $Id: vtkDICOMVolumeReader.h,v 1.2 2003/01/15 18:52:32 cpbotha Exp $
// class for reading off-line DICOM datasets

/* TODO
    - SliceThickness + SpacingBetweenSlices
  */

#ifndef vtkDICOMVolumeReader_h
#define vtkDICOMVolumeReader_h

#include <list>
#include <string>
#include <vector>
#include <vtkStructuredPointsSource.h>
#include "vtkdscasIOWin32Header.h"

// we need this else Visual C++ doesn't like our STL thingies
// in addition, this header should be included LAST in vtkDICOMVolumeReader.cxx
using namespace std;

//BTX
class DcmFileStream;
class DcmObject;
class DcmElement;
class DcmStack;

/**
 * Class used by vtkDICOMVolumeReader for keeping track of a bunch of dicom
 * images belonging to a dataset.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class dicom_file
{
public:
   string filename;
   DcmFileStream* filestream;
   DcmObject* fileformat;
   double SliceLocation;
   /**
    * With this operator defined, we can use stl to sort a vector of dicom_files
    * according to SliceLocation.  This is very handy.
    */
   bool operator < (const dicom_file& other_file) const
   {
      return SliceLocation < other_file.SliceLocation;
   }
};

class series_instance_misc_metadata_t
{
public:
   string StudyDescription;
   string ReferringPhysician;
};

/**
 * This class contains info about a specific series instance.  There will be one of
 * these for each instance scanned from the DICOM files.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class series_instance
{
public:
   string SeriesInstanceUID;
   double SliceThickness;
   double PixelSpacingx;
   double PixelSpacingy;
   unsigned short Rows;
   unsigned short Columns;
   unsigned short BitsAllocated;
   vector<dicom_file> dicom_files;
   series_instance_misc_metadata_t misc_metadata;
};
//ETX

/**
 * VTK class that makes use of DCMTK to enable dscas1 to read DICOM medical
 * volume data.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class VTK_DSCAS_IO_EXPORT vtkDICOMVolumeReader : public vtkStructuredPointsSource
{
protected:
   /// origin of the data, in true dimensions
   float DataOrigin[3];   
   /// Actual number of voxels on x, y and z axes.
   int DataDimensions[3];
   /// Millimetre spacing between voxels on x, y and z axes.
   float DataSpacing[3];
   /// List of all the DICOM image files that are going to be read by us.
   //BTX
   /// stores list of dicom filenames that are going to be read
   vector<string> dicom_filenames;
   /**
    * Filenames live here right after they've been add_dicom_filename()ed.  When
    * this->GetModified() is called, this buffer is compared with the internal dicom_filenames
    * vector.  If these are different, then MTime() is timestamped and the list is copied to
    * the internal vector.
    */
   vector<string> dicom_filenames_buffer;
   //ETX
   /// Index into internal vectors of filenames, each vector with unique SeriesInstanceUID
   int SeriesInstanceIdx;
   //BTX
   /** 
    * This list  contains series_instance classes, each with all the info pertaining to a
    * specific series instance.
    */
   list<series_instance> series_instances;
   //ETX
   double WindowCenter;
   double WindowWidth;

   vtkDICOMVolumeReader();
   ~vtkDICOMVolumeReader();

   void deinit_dcmtk(void);
   //BTX
   DcmElement* search_object(int group, int elem, DcmObject& haystack, DcmStack &stack);
   //ETX
   /**
     * Opens single DICOM file with filename in dfile.filename.
     * @param dfile dicom_file struct with valid filename field.  If successful, fileformat will
     * contain a valid DcmObject and filestream a valid DcmFileStream.
     * @return true if successful, false otherwise.
     */
   int OpenDCMFile(dicom_file& dfile);
   void ExecuteInformation(void);
   void ExecuteData(vtkDataObject* out);
   //BTX
   list<series_instance>::iterator find_si_iterator(int idx);
   //ETX

public:
   vtkTypeMacro(vtkDICOMVolumeReader, vtkStructuredPointsSource);
   static vtkDICOMVolumeReader* New();
   unsigned long GetMTime();
   /**
    * Add a single filename to the dicom_filenames_buffer list (i.e. not the
    * internal list).  The modified time is NOT set.  Whet GetMTime() is called,
    * the external and internal list will be checked for differences.  If they differ,
    * mtime will be timestamped.
    * @param dicom_filename Filename that you want to add to the buffer list
    * of dicom filenames.
    */
   void add_dicom_filename(const char* dicom_filename);
   /**
    * Remove all DICOM filenames from dicom_files_buffer and NOT the internal
    * list.  We also do not set the Modified time; this is so that if the user clears 
    * the list and re-adds the same filenames, the reader will still be up to date.  
    * This is (amongst other things) so that DSCAS3 doesn't let this reader re-read 
    * EVERY time that the user applies and the list of filenames is re-initialised...  
    * it also mimics the behaviour of a normal Set*() call in that if the parameter 
    * is the same is the internal version, nothing is done.
    */
   void clear_dicom_filenames(void);
   vtkGetMacro(WindowCenter, double);
   vtkGetMacro(WindowWidth, double);
   /// This will also make sure that ->Modified is set.
   vtkSetMacro(SeriesInstanceIdx, int);
   vtkGetMacro(SeriesInstanceIdx, int);
   vtkSetVector3Macro(DataOrigin,float);
   vtkGetVectorMacro(DataOrigin,float,3);
   vtkSetVector3Macro(DataDimensions,int);
   vtkGetVectorMacro(DataDimensions,int,3);
   vtkSetVector3Macro(DataSpacing,float);
   vtkGetVectorMacro(DataSpacing,float,3);

   const char *GetSeriesInstanceUID(void);
   const char *GetStudyDescription(void);
   const char *GetReferringPhysician(void);
};

#endif
