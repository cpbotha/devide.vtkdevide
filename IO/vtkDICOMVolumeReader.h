// vtkDICOMVolumeReader.h copyright (c) 2003 Charl P. Botha cpbotha@ieee.org
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id: vtkDICOMVolumeReader.h,v 1.15 2003/10/07 16:40:48 cpbotha Exp $
// class for reading off-line DICOM datasets

/*
 * This software is licensed exclusively for research use by Bart Kaptein
 * in the ModelBasedRSA package.  Any modifications made to this software
 * shall be sent to the author for possible inclusion in future versions.
 * Ownership and copyright of said modifications shall be ceded to the
 * author.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Example:
 * dreader = vtkDICOMVolumeReader::New();
 * dreader.add_dicom_filename('file1.dcm');
 * dreader.add_dicom_filename('file2.dcm');
 * ... rather do this with a loop: this code is robust enough to reject
 *     incorrect files (at read time) and also to discern automatically
 *     between different series
 * // this will automatically call an UpdateInformation()
 * maxsii = dreader.GetMaximumSeriesInstanceIdx();
 * for (int i; i <= maxsii; i++) {
 *   dreader.SetSeriesInstanceIdx(i);
 *   do_stuff_with(dreader.GetOutput());
 * }
 *
 * So, one uses SetSeriesInstanceIdx to select a different series from
 * the DICOM dataset.  The GetOutput() will "change" into a different
 * VTK dataset for each SeriesInstanceIdx.
 *
 * Because some DICOM datasets can be large, it is recommended that you
 * give progress feedback by making use of an Observer of the ProgressEvent.
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
//using namespace std;

//BTX
class DcmObject;
class DcmElement;
class DcmFileFormat;
class DcmStack;

/**
 * Class used by vtkDICOMVolumeReader for keeping track of a bunch of dicom
 * images belonging to a dataset.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class dicom_file
{
public:
   std::string filename;
   //DcmFileStream* filestream;
   DcmFileFormat* fileformat;
   double SliceLocation;
   /**
    * With this operator defined, we can use stl to sort a vector of
    * dicom_files according to SliceLocation.  This is very handy.
    */
   bool operator < (const dicom_file& other_file) const
   {
      return SliceLocation < other_file.SliceLocation;
   }
};

class series_instance_misc_metadata_t
{
public:
   std::string StudyDescription;
   std::string ReferringPhysician;
};

/**
 * This class contains info about a specific series instance.  There will be
 * one of these for each instance scanned from the DICOM files.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class series_instance
{
public:
   std::string SeriesInstanceUID;
   double SliceThickness;
   double SpacingBetweenSlices;
   double PixelSpacingx;
   double PixelSpacingy;
   unsigned short Rows;
   unsigned short Columns;
   unsigned short BitsAllocated;
   unsigned short PixelRepresentation;
   unsigned short BitsStored;
   std::vector<dicom_file> dicom_files;
   series_instance_misc_metadata_t misc_metadata;
};
//ETX

/**
 * VTK class that makes use of DCMTK to enable dscas1 to read DICOM medical
 * volume data.
 * @author Charl P. Botha <cpbotha@ieee.org>
 */
class VTK_DSCAS_IO_EXPORT vtkDICOMVolumeReader :
public vtkStructuredPointsSource
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
   std::vector<std::string> dicom_filenames;
   /**
    * Filenames live here right after they've been add_dicom_filename()ed.
    * When this->GetModified() is called, this buffer is compared with the
    * internal dicom_filenames vector.  If these are different, then MTime()
    * is timestamped and the list is copied to the internal vector.
    */
   std::vector<std::string> dicom_filenames_buffer;
   //ETX
   /// Index into internal vectors of filenames, each vector with unique
   /// SeriesInstanceUID
   int SeriesInstanceIdx;
   //BTX
   /** 
    * This list  contains series_instance classes, each with all the info
    * pertaining to a specific series instance.
    */
   std::list<series_instance> series_instances;
   //ETX
   double WindowCenter;
   double WindowWidth;

   /// When on, this reader will be far more lenient with the files it
   /// has to parse.
   int Leniency;

   vtkDICOMVolumeReader();
   ~vtkDICOMVolumeReader();

   void deinit_dcmtk(void);
   //BTX
   DcmElement* search_object(int group, int elem, DcmObject& haystack, DcmStack &stack);
   //ETX
   /**
     * Opens single DICOM file with filename in dfile.filename.
     * @param dfile dicom_file struct with valid filename field.  If
     * successful, fileformat will contain a valid DcmObject and filestream
     * a valid DcmFileStream.
     * @return true if successful, false otherwise.
     */
   int OpenDCMFile(dicom_file& dfile);
   void ExecuteInformation(void);
   void ExecuteData(vtkDataObject* out);
   //BTX
   std::list<series_instance>::iterator find_si_iterator(int idx);
   //ETX

public:
   vtkTypeMacro(vtkDICOMVolumeReader, vtkStructuredPointsSource);
   static vtkDICOMVolumeReader* New();
   unsigned long GetMTime();
   /**
    * Add a single filename to the dicom_filenames_buffer list (i.e. not the
    * internal list).  The modified time is NOT set.  Whet GetMTime() is
    * called, the external and internal list will be checked for differences.
    * If they differ, mtime will be timestamped.
    * @param dicom_filename Filename that you want to add to the buffer list
    * of dicom filenames.
    */
   void add_dicom_filename(const char* dicom_filename);
   /**
    * Remove all DICOM filenames from dicom_files_buffer and NOT the internal
    * list.  We also do not set the Modified time; this is so that if the user
    * clears the list and re-adds the same filenames, the reader will still be
    * up to date.  This is (amongst other things) so that DSCAS3 doesn't let
    * this reader re-read EVERY time that the user applies and the list of
    * filenames is re-initialised...  it also mimics the behaviour of a normal
    * Set*() call in that if the parameter is the same is the internal version,
    * nothing is done.
    */
   void clear_dicom_filenames(void);
   /**
    * Return number of filenames in dicom_filenames_buffer list.
    */
   int get_number_of_dicom_filenames(void);
   /**
    * Return idx'th filename in dicom_filenames_buffer list.  If idx is
    * invalid, returns NULL.
    */
   const char *get_dicom_filename(int idx);
   vtkGetMacro(WindowCenter, double);
   vtkGetMacro(WindowWidth, double);

   vtkSetMacro(Leniency, int);
   vtkGetMacro(Leniency, int);
   vtkBooleanMacro(Leniency, int);

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

   /**
    * Return the maximum SeriesInstanceIdx.
    */
   int GetMaximumSeriesInstanceIdx(void);
};

#endif
