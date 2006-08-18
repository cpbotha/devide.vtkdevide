// vtkDICOMVolumeReader.cxx copyright (c) 2001-2006 Charl P. Botha
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id$
// class for reading off-line DICOM datasets

/*
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

/* Keep logic as it is now, but use DicomImage to read PixelData
   via the DicomImage(DcmObject...) form.

   See: http://www.archivaty.com/forums/about68765.html
        http://support.dcmtk.org/docs/file_faq.html (compression)

   Maybe not, only in DCMTK 3.5.4 does DicomImage give access to
   the underlying pixeldata.  Let's wait this out.
*/

#if !defined(WIN32)
#define HAVE_CONFIG_H
#endif
#include <osconfig.h>
#include <dcmdata/dctk.h>

//#include <dcmjpeg/djdecode.h>
//#include <dcmdata/dcrledrg.h>

#include <algorithm>
#include <math.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkStructuredPoints.h>

#include "vtkDICOMVolumeReader.h"

#define INTROUND(a) (int)(a + 0.5f)

// ##########################################################################
vtkDICOMVolumeReader::vtkDICOMVolumeReader()
{
  // try this again for jpeg stuff
  //DJDecoderRegistration::registerCodecs();
  
  // this should be a sane default
  this->SeriesInstanceIdx = 0;
  this->series_instances.clear();
  this->dicom_filenames.clear();
  this->dicom_filenames_buffer.clear();

  // we have no valid spacing yet, so we show it this way
  DataSpacing[0] = DataSpacing[1] = DataSpacing[2] = 0.0;
  DataOrigin[0] = DataOrigin[1] = DataOrigin[2] = 0.0;
  DataDimensions[0] = DataDimensions[1] = DataDimensions[2] = 0;
  this->WindowWidth = this->WindowCenter = 0.0;

  this->LeniencyOff();
  this->EstimateSliceThicknessOn();
}

// ##########################################################################
vtkDICOMVolumeReader::~vtkDICOMVolumeReader()
{
  deinit_dcmtk();
}

// ##########################################################################
vtkDICOMVolumeReader* vtkDICOMVolumeReader::New()
{
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkDICOMVolumeReader");
  if (ret)
    {
    return(vtkDICOMVolumeReader*)ret;
    }
  return new vtkDICOMVolumeReader;
}

// ##########################################################################
void vtkDICOMVolumeReader::deinit_dcmtk(void)
{
  vtkDebugMacro(<<"vtkDICOMVolumeReader::deinit_dcmtk() - START.");
  std::list<series_instance>::iterator si_iterator;
  // go through all series instances
  for (si_iterator = series_instances.begin(); si_iterator != series_instances.end(); si_iterator++)
    {
    // in each series instance, call dtor of fileformats
    std::vector<dicom_file>* dicom_files_p = &((*si_iterator).dicom_files);
    for (unsigned i = 0; i < dicom_files_p->size(); i++)
      {
      if ((*dicom_files_p)[i].fileformat) delete (*dicom_files_p)[i].fileformat;
      }
    // clear the vectors to play it safe
    dicom_files_p->clear();
    }
   
  // then empty the list itself
  series_instances.clear();

  // deinit codecs
  //DJDecoderRegistration::cleanup();

  // we are very seriously modified
  this->Modified();
  vtkDebugMacro(<<"vtkDICOMVolumeReader::deinit_dcmtk() - DONE.");
}

// ##########################################################################
DcmElement* vtkDICOMVolumeReader::search_object(int group, int elem, DcmObject& haystack, DcmStack &stack)
{
  DcmTagKey search_key;
  search_key.set(group,elem);
  if (haystack.search(search_key, stack, ESM_fromHere, OFTrue) != EC_Normal)
    return NULL;


  else
    {
    // we only want elements at top level (i.e. not nested) -- for elements at top level
    // that have been found, the stack is 2 deep (one item for the list itself, one
    // for the found element)
    if (stack.card() > 2)
      {
      E_Condition sres;
      // so we continue searching (using ESM_afterStackTop) until we've found the top
      // level element that we were looking for
      do
        {
        sres = haystack.search(search_key, stack, ESM_afterStackTop, OFTrue);
        } while (sres == EC_Normal && stack.card() > 2);
      if (sres != EC_Normal || stack.card() > 2)
        return NULL;
      }
    // if we haven't returned by now, it means we've found our thing
    return(DcmElement*)stack.top();
    }
}

// ##########################################################################
int vtkDICOMVolumeReader::OpenDCMFile(dicom_file& dfile)
{
  vtkDebugMacro(<<"vtkDICOMVolumeReader::OpenDCMFile(" << dfile.filename.c_str() << ") - START.");

  // we have the stream, now we instantiate a DICOM fileformat abstraction to read and interpret the data
  // DcmFileFormat is a DcmSequence (or something) containing a list of DcmItems, each representing a separate
  // (group,element) data tag
  dfile.fileformat = new DcmFileFormat();
  E_TransferSyntax xfer = EXS_Unknown;
  OFCondition cond = dfile.fileformat->loadFile(dfile.filename.c_str());
    //dfile.filename.c_str(), xfer, EGL_noChange, DCM_MaxReadLength, false);
  if (! cond.good())
    {
    vtkErrorMacro(<<"Unable to open/read DICOM file " << dfile.filename.c_str());
    return 0;
    }

  vtkDebugMacro(<<"vtkDICOMVolumeReader::OpenDCMFile() - STOP.");
  return 1;
}

// ##########################################################################
void vtkDICOMVolumeReader::ExecuteInformation(void)
{
  vtkDebugMacro(<<"vtkDICOMVolumeReader::ExecuteInformation() - START");

  // go through whole list of files, get SeriesInstanceUIDs, separate
  // files into different lists

  // we have to make sure the DICOM dictionary is loaded (this
  // requires a functioning DICOM installation)
  if (!dcmDataDict.isDictionaryLoaded())
    {
    vtkErrorMacro("ExecuteInformation(): The DICOM dictionary isn't loaded!");
    return;
    }

  // check that the user/programmer has given us (effectively) a list
  // of dicom image files
  if (dicom_filenames.size() <= 0)
    {
    return;
    }

  // make sure we start from scratch
  this->deinit_dcmtk();

  
  char* SeriesInstanceUID_cp;
  std::list<series_instance>::iterator si_iterator;

  for (unsigned i = 0; i < dicom_filenames.size(); i++)
    {
    dicom_file temp_dicom_file;
    temp_dicom_file.fileformat = NULL;
      
    temp_dicom_file.filename = dicom_filenames[i];
    if (!this->OpenDCMFile(temp_dicom_file))
      {
      // if we can't read one of the DICOM files, just continue to the next one
      if (temp_dicom_file.fileformat)
        delete temp_dicom_file.fileformat;
      continue;
      }

    // get SeriesInstanceUID
    // ---------------------------------------------------------------------
    SeriesInstanceUID_cp = NULL;
    DcmStack SeriesInstanceStack;
    DcmElement* SeriesInstanceUID_obj = search_object(0x0020,0x000e, *(temp_dicom_file.fileformat), SeriesInstanceStack); // SeriesInstanceUID
    if (!SeriesInstanceUID_obj || (SeriesInstanceUID_obj->getString(SeriesInstanceUID_cp) != EC_Normal) || !SeriesInstanceUID_cp)
      {
      vtkWarningMacro(<< "vtkDICOMVolumeReader::ExecuteInformation() - Unable to read SeriesInstanceUID from " << temp_dicom_file.filename.c_str() << ".");
      // we can just skip this file, starting at the next iteration
      // first make sure what's there now gets deallocation/de-initialised
      if (temp_dicom_file.fileformat)
        delete temp_dicom_file.fileformat;
      // is this frowned upon like goto?
      continue;
      }

    std::string SeriesInstanceUID_str = std::string(SeriesInstanceUID_cp);

    // check if this series_instance exists...
    bool si_found = false;
    for (si_iterator = series_instances.begin();
	 !si_found && si_iterator != series_instances.end(); si_iterator++)
      {
      si_found = (*si_iterator).SeriesInstanceUID == SeriesInstanceUID_str;
      }

    if (si_found)
      {
      // the for loop always gets in one increment at the end of the loop!
      si_iterator--;
      } // if (si_found)
    else
      {
      vtkDebugMacro(<< "New SeriesInstance "
		    << SeriesInstanceUID_str.c_str() << " found.");
      
      // do everything we need to do with a new instance (patient
      // info, dataset info, ImageOrientationPatent, etc.)
      // for now we assume that everything within an instance has to
      // be the same
      series_instance temp_series_instance;
      temp_series_instance.SeriesInstanceUID = SeriesInstanceUID_str;

      // now make sure SliceThickness, PixelSpacing, Rows, Columns are
      // all present
      double temp_SliceThickness, temp_SpacingBetweenSlices;
      double temp_PixelSpacingx, temp_PixelSpacingy;
      unsigned short temp_Rows, temp_Columns;
      unsigned short temp_BitsAllocated;
      unsigned short temp_PixelRepresentation;

      DcmStack SliceThickness_stack;

      DcmElement* SliceThickness_obj = search_object(
	0x0018, 0x0050, *(temp_dicom_file.fileformat), SliceThickness_stack); 
      if (!SliceThickness_obj ||
	  SliceThickness_obj->getFloat64(temp_SliceThickness) != EC_Normal)
        {
        if (!this->GetLeniency())
          {
          vtkErrorMacro(
	    <<"::ExecuteInfo() - could not read SliceThickness from "
	    << temp_dicom_file.filename.c_str() << ", ignoring file.");
          if (temp_dicom_file.fileformat)
            delete temp_dicom_file.fileformat;
          continue;
          }
        else
          {
          vtkWarningMacro(
	    <<"::ExecuteInfo() - could not read SliceThickness from "
	    << temp_dicom_file.filename.c_str() << ", assuming 1mm.");
          temp_SliceThickness = 1;
          }
        }

      DcmStack SpacingBetweenSlices_stack;
      DcmElement* SpacingBetweenSlices_obj = search_object(
	0x0018, 0x0088, *(temp_dicom_file.fileformat),
	SpacingBetweenSlices_stack); // SpacingBetweenSlices
      
      if (!SpacingBetweenSlices_obj ||
	  SpacingBetweenSlices_obj->getFloat64(
	    temp_SpacingBetweenSlices) != EC_Normal)
        {
        // this is normal: only some MRI datasets make use of this
        temp_SpacingBetweenSlices = -1;
        }
      else
        {
        // sometimes we find a negative SpacingBetweenSlices (see Rommens)
        // in this case, we do like Atamai does and fabs it.
        temp_SpacingBetweenSlices = fabs(temp_SpacingBetweenSlices);
        }

      // ---------------------------------------------------------------------
      // we read the ImageOrientationPatient as we want to do
      // something like:
      // http://www.itk.org/pipermail/insight-users/2003-September/004762.html
      DcmStack ImageOrientationPatient_stack;
      DcmElement* ImageOrientationPatient_obj = search_object(
	0x0020, 0x0037, *(temp_dicom_file.fileformat),
	ImageOrientationPatient_stack);

      double temp_iop[6];


      // ARGH! ImageOrientationPatient_obj->getFloat64Array() errors.
      // we HAVE to do it one by one with getFloat64Array(var, idx)
      bool iop_success = true;
      if (!ImageOrientationPatient_obj)
	{
	iop_success = false;
	}
      else
	{
	for (int i = 0; i < 6  && iop_success; i++)
	  {
	  if (ImageOrientationPatient_obj->getFloat64(
		temp_iop[i], i) != EC_Normal)
	    {
	    iop_success = false;
	    }
	  }
	}

      if (!iop_success)
      {

	vtkErrorMacro(
	  <<"::ExecuteInfo() - could not read ImageOrientationPatient from "
	  << temp_dicom_file.filename.c_str() << ", ignoring file.");

	if (temp_dicom_file.fileformat)
	  delete temp_dicom_file.fileformat;
	continue;
      }


      DcmStack PixelSpacing_stack;
      DcmElement* PixelSpacing_obj = search_object(
	0x0028, 0x0030, *(temp_dicom_file.fileformat), PixelSpacing_stack);
      
      if (!PixelSpacing_obj || 
          PixelSpacing_obj->getFloat64(temp_PixelSpacingx,0) != EC_Normal ||
          PixelSpacing_obj->getFloat64(temp_PixelSpacingy,1) != EC_Normal)
        {
        if (!this->GetLeniency())
          {
          vtkErrorMacro(<<"::ExecuteInfo() - could not read PixelSpacing from " << temp_dicom_file.filename.c_str() << ", ignoring file.");
          if (temp_dicom_file.fileformat)
            delete temp_dicom_file.fileformat;
          continue;
          }
        else
          {
          vtkWarningMacro(<<"::ExecuteInfo() - could not read PixelSpacing from " << temp_dicom_file.filename.c_str() << ", assuming 1mm.");
          temp_PixelSpacingx = 1;
          temp_PixelSpacingy = 1;
          }
        }

      DcmStack Rows_stack, Columns_stack;
      DcmElement* Rows_obj = search_object(0x0028, 0x0010, *(temp_dicom_file.fileformat), Rows_stack); // Rows
      DcmElement* Columns_obj = search_object(0x0028, 0x0011, *(temp_dicom_file.fileformat), Columns_stack); // Columns
      if (!Rows_obj || !Columns_obj || Rows_obj->getUint16(temp_Rows) != EC_Normal || Columns_obj->getUint16(temp_Columns) != EC_Normal)
        {
        vtkErrorMacro(<<"::ExecuteInfo() - could not read Rows/Colums from " << temp_dicom_file.filename.c_str() << ", ignoring file.");
        if (temp_dicom_file.fileformat)
          delete temp_dicom_file.fileformat;
        continue;
        }

      DcmStack BitsAllocated_stack;
      DcmElement* BitsAllocated_obj = search_object(0x0028, 0x0100, *(temp_dicom_file.fileformat), BitsAllocated_stack); // BitsAllocated
      if (!BitsAllocated_obj || BitsAllocated_obj->getUint16(temp_BitsAllocated) != EC_Normal)
        {
        vtkWarningMacro(<<"Could not read BitsAllocated from " << temp_dicom_file.filename.c_str() << ", assuming 16.");
        temp_BitsAllocated = 16;
        }

      // we need this so we know if our output has to be short or unsigned short
      DcmStack PixelRepresentation_stack;
      DcmElement* PixelRepresentation_obj = search_object(0x0028, 0x0103, *(temp_dicom_file.fileformat), PixelRepresentation_stack); // PixelRepresentation
      if (!PixelRepresentation_obj || PixelRepresentation_obj->getUint16(temp_PixelRepresentation) != EC_Normal)
        {
        vtkWarningMacro(<<"Could not read PixelRepresentation from " << temp_dicom_file.filename.c_str() << ", assuming 1 (signed).");
        temp_PixelRepresentation = 1;
        }

      DcmStack BitsStored_stack;
      unsigned short temp_BitsStored;
      DcmElement* BitsStored_obj = search_object(0x0028, 0x0101, *(temp_dicom_file.fileformat), BitsStored_stack);
      if (!BitsStored_obj || BitsStored_obj->getUint16(temp_BitsStored) != EC_Normal)
        {
        vtkWarningMacro(<<"Could not read BitsStored from " << temp_dicom_file.filename.c_str() << ", assuming 12 (normal for CT).");
        temp_BitsStored = 12;
        }


      // then some arb shit...
      char *StudyDescription_cp;
      DcmStack StudyDescription_stack;
      DcmElement *StudyDescription_obj = search_object(0x0008,0x1030, *(temp_dicom_file.fileformat), StudyDescription_stack); // StudyDescription
      if (!StudyDescription_obj || StudyDescription_obj->getString(StudyDescription_cp) != EC_Normal || StudyDescription_cp == NULL)
        temp_series_instance.misc_metadata.StudyDescription = std::string("N/A");
      else
        temp_series_instance.misc_metadata.StudyDescription = std::string(StudyDescription_cp);

      char *ReferringPhysician_cp;
      DcmStack ReferringPhysician_stack;
      DcmElement *ReferringPhysician_obj = search_object(0x0008,0x0090, *(temp_dicom_file.fileformat), ReferringPhysician_stack); // ReferringPhysician
      if (!ReferringPhysician_obj || ReferringPhysician_obj->getString(ReferringPhysician_cp) != EC_Normal || ReferringPhysician_cp == NULL)
        temp_series_instance.misc_metadata.ReferringPhysician = std::string("N/A");
      else
        temp_series_instance.misc_metadata.ReferringPhysician = std::string(ReferringPhysician_cp);


      // store all the other thingies we've just extracted
      temp_series_instance.SliceThickness = temp_SliceThickness;
      temp_series_instance.SpacingBetweenSlices = temp_SpacingBetweenSlices;

      for (int i = 0; i < 6; i++)
	temp_series_instance.ImageOrientationPatient[i] = temp_iop[i];

      // calculate iop normal (z-normal) and store it
      temp_series_instance.IOPNormal[0] = temp_iop[1]*temp_iop[5] - \
	temp_iop[2]*temp_iop[4];
      temp_series_instance.IOPNormal[1] = temp_iop[2]*temp_iop[3] - \
	temp_iop[0]*temp_iop[5];
      temp_series_instance.IOPNormal[2] = temp_iop[0]*temp_iop[4] - \
	temp_iop[1]*temp_iop[3];
      
      temp_series_instance.PixelSpacingx = temp_PixelSpacingx;
      temp_series_instance.PixelSpacingy = temp_PixelSpacingy;
      temp_series_instance.Rows = temp_Rows;
      temp_series_instance.Columns = temp_Columns;
      temp_series_instance.BitsAllocated = temp_BitsAllocated;
      temp_series_instance.PixelRepresentation = temp_PixelRepresentation;
      temp_series_instance.BitsStored = temp_BitsStored;

      // make sure the dicom_files vector is empty
      temp_series_instance.dicom_files.clear();
      // we now should have a complete series_instance, so we add it
      // to the linked list
      series_instances.push_back(temp_series_instance);
      // we set si_iterator to the just push_back()ed series_instance
      // so the following code doesn't have to be conditional
      si_iterator = series_instances.end();
      si_iterator--;
      } // if (si_found) ... else (new SI, iow)

    // we get the found object (DcmObject -> DcmElement -> ?, in this
    // case probably DcmFloatingPointDouble) from the stack


    DcmStack ImagePositionPatient_stack;
    DcmElement* ImagePositionPatient_obj = search_object(
      0x0020, 0x0032, *(temp_dicom_file.fileformat),
      ImagePositionPatient_stack);

    double temp_ipp[3];
    double distance = 0;

    bool ipp_success = true;
    
    if (!ImagePositionPatient_obj)
      {
      ipp_success = false;
      }

    else
      {
      for (int i = 0; i < 3 && ipp_success; i++)
	{
	if (ImagePositionPatient_obj->getFloat64(temp_ipp[i],i) != EC_Normal)
	  {
	  ipp_success = false;
	  }
	}
      }

    if (!ipp_success)
      {
      vtkWarningMacro(
	<<"::ExecuteInfo() - could not read ImagePositionPatient from "
	<< temp_dicom_file.filename.c_str() << ", using SliceLocation.");

      double temp_SliceLocation;
      DcmStack SliceLocation_stack;
      DcmElement* SliceLocation_obj = search_object(
	0x0020, 0x1041, *(temp_dicom_file.fileformat), SliceLocation_stack);
    
      if (!SliceLocation_obj ||
	  (SliceLocation_obj->getFloat64(temp_SliceLocation) != EC_Normal))
	{
	vtkWarningMacro(
	  <<"::ExecuteInfo() - could not read SliceLocation from "
	  << temp_dicom_file.filename.c_str() << ", using 0.");
	temp_SliceLocation = 0;
	}

      distance = temp_SliceLocation;

      }
    else
      {
      // this is the Right Way(tm) of calcing the distance!
      for (int i = 0; i < 3; i++)
	distance += (*si_iterator).IOPNormal[i] * temp_ipp[i];

      }

    temp_dicom_file.distance = distance;
    
      
    // si_iterator is the series_instance where the current
    // SeriesInstanceUID was found, so we need to store the dicom_file
    // object
    (*si_iterator).dicom_files.push_back(temp_dicom_file);

    } // for (i = 0; i < dicom_filenames.size(); i++) ...

  for (si_iterator = series_instances.begin();
       si_iterator != series_instances.end(); si_iterator++)
    {
    // sort dicom_files according to "distance", hmmmmkay?
    std::sort((*si_iterator).dicom_files.begin(),
	      (*si_iterator).dicom_files.end());
    
    vtkDebugMacro(
      <<"Sorted SeriesInstanceUID "
      << (*si_iterator).SeriesInstanceUID.c_str() << ", "
      << (*si_iterator).BitsAllocated << ", "
      << (*si_iterator).dicom_files.size() << " files.");

    // also find "most popular" slice thickness according to "distance"
    std::vector<dicom_file>* dicom_files_p = &((*si_iterator).dicom_files);
    std::vector<double > diffs;  diffs.clear();
    std::vector<int > diffCounts;  diffCounts.clear();
    bool diffFound;
    double tempDiff;
    for (int sliceIdx = 1; sliceIdx < dicom_files_p->size(); sliceIdx++)
      {
      tempDiff =
        (*dicom_files_p)[sliceIdx].distance -
        (*dicom_files_p)[sliceIdx-1].distance;

      // now find tempDiff in diffs
      diffFound = false;
      int i;
      for (i = 0; i < diffs.size() && !diffFound; i++)
        {
        if (fabs(diffs[i] - tempDiff) < 0.01)
          {
          diffFound = true;
          }
        }

      if (diffFound)
        {
        i--;
        diffCounts[i]++;
        }
      else
        {
        diffs.push_back(tempDiff);
        diffCounts.push_back(1);
        }
      }

    // now go through diffs and diffCounts, find difference with
    // largest number of occurrences
    int maxi = -1;
    int maxCount = 0;
    for (int i = 0; i < diffCounts.size(); i++)
      {
      if (diffs[i] > maxCount)
        {
        maxi = i;
        maxCount = diffCounts[i];
        }
      }

    if (maxi >= 0)
      {
      (*si_iterator).EstimatedThickness = diffs[maxi];
      }
    else
      {
      // this should never happen, but we have a contingency plan
      (*si_iterator).EstimatedThickness = (*si_iterator).SliceThickness;
      }
    }

  // ------------------------------------------------------------------------

  // the following depends on the currently selected SeriesInstanceIdx
  si_iterator = this->find_si_iterator(this->SeriesInstanceIdx);
  // if it's not found, si_iterator == series_instances.end()
  // we have to yield something, so let's show the last series instance then
  if (si_iterator == series_instances.end())
    {
    si_iterator--;
    }

  // after this loop, we have reached either the SeriesInstanceIdx'th
  // series_instance or the highest, whichever comes first and we can
  // set the data that has to be valid for the current series

  // transfer IOP to ivar
  for (int i = 0; i < 6; i++)
    {
    this->ImageOrientationPatient[i] =
      (*si_iterator).ImageOrientationPatient[i];
    }
  
  this->DataDimensions[0] = (*si_iterator).Columns;
  this->DataDimensions[1] = (*si_iterator).Rows;
  this->DataDimensions[2] = (*si_iterator).dicom_files.size();
  this->DataSpacing[0] = (*si_iterator).PixelSpacingx;
  this->DataSpacing[1] = (*si_iterator).PixelSpacingy;
   
  // SpacingBetweenSlices is preferred for the axial spacing.  It's
  // mostly only present in MRI datasets.  CT datasets have
  // SliceThickness, but it seems that SliceThickness often does not
  // concur with IPP-derived distance.  So, if this is the case, we
  // use the IPP/distance derived thickness (i.e. EstimatedThickness)
  if ((*si_iterator).SpacingBetweenSlices > 0)
    {
    DataSpacing[2] = (*si_iterator).SpacingBetweenSlices;
    }
  else
    {
    if (this->EstimateSliceThickness)
      {
      DataSpacing[2] = (*si_iterator).EstimatedThickness;
      }
    else
      {
      DataSpacing[2] = (*si_iterator).SliceThickness;
      }

    }

  // get pointer to our output vtkStructuredPoints
  vtkStructuredPoints* output = this->GetOutput();

  // now tell it how many voxels it's going to have in xyz directions
  // confusingly, this call actually sets up the extent (and when getdimensions
  // are called, the result is also calculated from the extent)
  output->SetDimensions(DataDimensions);
  // and tell it the millimetre spacing between voxels
  output->SetSpacing(DataSpacing[0], DataSpacing[1], DataSpacing[2]);
  // and the origin (in physical coordinates)
  output->SetOrigin(DataOrigin[0], DataOrigin[1], DataOrigin[2]);
  // and setup the WholeExtent to be the same as the extent (Extent is
  // what's in mem, WholeExtent is everything)
  output->SetWholeExtent(
    0, DataDimensions[0]-1, 0, DataDimensions[1]-1, 0, DataDimensions[2]-1);

  if ((*si_iterator).PixelRepresentation == 0 &&
      (*si_iterator).BitsStored == 16)
    {
    // this is a special case... by definition, it has to be the full
    // 16 bits, so we have to use unsigneds.  You better hope that the
    // rescaleIntercept is not negative, then this will break!
    output->SetScalarType(VTK_UNSIGNED_SHORT);
    }
  else
    {
    // we try to have signed data wherever possible.  Even if the
    // pixeldata is unsigned, the RescaleIntercept (or slope) could
    // make us negative!
    output->SetScalarType(VTK_SHORT);
    }

  //output->SetScalarType((*si_iterator).BitsAllocated == 8 ?
  //VTK_UNSIGNED_CHAR : VTK_SHORT); and we will have only one scalar
  //component
  output->SetNumberOfScalarComponents(1);
  
}

// ##########################################################################
void vtkDICOMVolumeReader::ExecuteData(vtkDataObject* out)
{
  vtkDebugMacro(<<"vtkDICOMVolumeReader::ExecuteData() - START");

  // we have to make sure the DICOM dictionary is loaded (this
  // requires a functioning DICOM installation)
  if (!dcmDataDict.isDictionaryLoaded())
    {
    vtkErrorMacro(<<"vtkDICOMVolumeReader::ExecuteData() - the DICOM dictionary isn't loaded!");
    return;
    }

  // check that the user/programmer has given us (effectively) a list
  // of dicom image files
  if (dicom_filenames.size() <= 0)
    {
    return;
    }

  vtkStructuredPoints* output = vtkStructuredPoints::SafeDownCast(out);
  if (!output)
    {
    vtkErrorMacro(<<"vtkDICOMVolumeReader::ExecuteData with non-vtkStructuredPoints output!");
    return;
    }

  // we're going to read the whole dataset, so make sure the extent is
  // good we could just do SetExtent(output->GetUpdateExtent()), but
  // then we would have to be more clever with our actual reading
  // routines, i.e. only read the data necessary to create the
  // UpdateExtent in the output
  output->SetExtent(output->GetWholeExtent());
  int *e = output->GetWholeExtent();
  cout << e[0] << ":" << e[1] << ":"
       << e[2] << ":" << e[3] << ":"
       << e[4] << ":" << e[5] << endl;
  // AllocateScalars will make use of the information set on the output by 
  // ExecuteInformation() to perform the necessary allocation
  output->AllocateScalars();
  cout << "MAXID: " << output->GetPointData()->GetScalars()->GetMaxId()
       << endl;


  // find a pointer that we can use
  short* spixels = (short*)output->GetScalarPointer();
  unsigned short *uspixels = (unsigned short*)output->GetScalarPointer();
   
  int numpixels = DataDimensions[0] * DataDimensions[1] * DataDimensions[2];
  int numpixels_perslice = DataDimensions[0] * DataDimensions[1];

  // find the iterator corresponding to the current seriesinstance
  std::list<series_instance>::iterator si_iterator = this->find_si_iterator(
    this->SeriesInstanceIdx);

  if (si_iterator == series_instances.end())
    si_iterator--;

  std::vector<dicom_file>* dicom_files_p = &((*si_iterator).dicom_files);
   
  // it's all sorted, so we can just load the pixel data and stuff it
  // a slice at a time into our allocated memory
  int last_progress_i = -1;
  for (unsigned i = 0; i < dicom_files_p->size(); i++)
    {

    double progress = (double)i / (double)dicom_files_p->size();
    // stupid way of only updating every 10 percent
    if ((progress / 0.1 - floor(progress / 0.1)) < 0.1 && (int)i != last_progress_i)
      {
      this->UpdateProgress(progress);
      last_progress_i = i;
      }

    DcmStack PixelData_stack;
    DcmElement* PixelData_obj = search_object(0x7fe0, 0x0010, *((*dicom_files_p)[i].fileformat), PixelData_stack);
    if (!PixelData_obj)
      {
      vtkErrorMacro(<<"vtkDICOMVolumeReader::ExecuteData() - could not extract PixelData (7fe0,0010)!");
      return;
      }

    // make sure it's OW (for now we only support OW)
    if (PixelData_obj->getTag().getEVR() != EVR_OW) // get Value Representation Enumeration
      {
      vtkErrorMacro(<<"vtkDICOMVolumeReader::ExecuteData() - PixelData not of type OW.");
      return;
      }

    // ------------------------------------------------------------------
    // Now we have to figure out the pixel representation and packing
    // so we know how to get out the SV (stored values)
    // ------------------------------------------------------------------
    DcmStack BitsStored_stack, HighBit_stack, PixelRepresentation_stack;
    DcmElement* BitsStored_obj = search_object(0x0028, 0x0101, *((*dicom_files_p)[i].fileformat), BitsStored_stack);
    DcmElement* HighBit_obj = search_object(0x0028, 0x0102, *((*dicom_files_p)[i].fileformat), HighBit_stack);
    DcmElement* PixelRepresentation_obj = search_object(0x0028, 0x0103, *((*dicom_files_p)[i].fileformat), PixelRepresentation_stack); // PixelRepresentation
    unsigned short temp_BitsStored, temp_HighBit, temp_PixelRepresentation;
    if (!BitsStored_obj || !HighBit_obj || !PixelRepresentation_obj || 
        BitsStored_obj->getUint16(temp_BitsStored) != EC_Normal ||
        HighBit_obj->getUint16(temp_HighBit) != EC_Normal ||
        PixelRepresentation_obj->getUint16(temp_PixelRepresentation) != EC_Normal)
      {
      vtkErrorMacro(<<"vtkDICOMVolumeReader::Excecute() - Could not get pixel representation parameters.");
      return;
      }
    unsigned rshift_factor = temp_HighBit + 1 - temp_BitsStored;
    unsigned validbits_mask = 0, highbit_mask = 1;
    for (int bit_i = 0; bit_i < temp_BitsStored; bit_i++)
      {
      highbit_mask <<= 1;
      validbits_mask <<= 1;
      // this will end up with a bit for each valid bit
      validbits_mask |= 1; // bytefactorl
      }
    // we want this bit right on the high bit for two's complement
    highbit_mask  >>= 1;
    //unsigned validbitsbuthigh_mask = validbits_mask >> 1; // bytemask

    // ------------------------------------------------------------------
    // And then we extract the Rescale intercept and slope, because we
    // need these to get the Hounsfield Units from extracted SVs
    // (stored values): HU = slope * SV + intercept
    // ------------------------------------------------------------------
    
    DcmStack RescaleIntercept_stack, RescaleSlope_stack;
    DcmElement* RescaleIntercept_obj = search_object(0x0028, 0x1052, *((*dicom_files_p)[i].fileformat), RescaleIntercept_stack); // RescaleIntercept
    DcmElement* RescaleSlope_obj = search_object(0x0028, 0x1053, *((*dicom_files_p)[i].fileformat), RescaleSlope_stack); // RescaleIntercept
    double temp_RescaleIntercept, temp_RescaleSlope;
    if (!RescaleIntercept_obj || !RescaleSlope_obj || 
        RescaleIntercept_obj->getFloat64(temp_RescaleIntercept) != EC_Normal ||
        RescaleSlope_obj->getFloat64(temp_RescaleSlope) != EC_Normal)
      {
      // it's possible that the DICOM does not have these two tags (e.g. some MRIs) in which case
      // HU is not pertinent in anycase; so we fake these constants
      temp_RescaleSlope = 1.0;
      temp_RescaleIntercept = 0.0;
      }

    // ---------------------------------------------------------------------
    // We extract WindowWidth and WindowCenter as well.  We don't need
    // these for getting the HU values, but our vtkDEVIDEVolumeReader
    // heritage requires that we have these available.
    // ---------------------------------------------------------------------
    
    DcmStack WindowCenter_stack, WindowWidth_stack;
    DcmElement* WindowCenter_obj = search_object(
      0x0028, 0x1050, *((*dicom_files_p)[i].fileformat), WindowCenter_stack);
    DcmElement* WindowWidth_obj = search_object(
      0x0028, 0x1051, *((*dicom_files_p)[i].fileformat), WindowWidth_stack);
    
    if (!WindowCenter_obj || !WindowWidth_obj || 
        WindowCenter_obj->getFloat64(this->WindowCenter) != EC_Normal ||
        WindowWidth_obj->getFloat64(this->WindowWidth) != EC_Normal)
      {
      vtkErrorMacro(<<"Could not extract Window{Center,Width}.");
      return;
      }

    //unsigned short* dicom_pixeldata_pointer;      
    unsigned short *dicom_pixeldata_pointer_us;
    unsigned char *dicom_pixeldata_pointer_uc;
    PixelData_obj->getUint16Array(dicom_pixeldata_pointer_us);
    PixelData_obj->getUint8Array(dicom_pixeldata_pointer_uc);

    // I need to document this bit of logic, because I don't get it anymore
    if ((*si_iterator).BitsAllocated == 8 && !dicom_pixeldata_pointer_us ||
        (*si_iterator).BitsAllocated != 8 &&  !dicom_pixeldata_pointer_uc)
      {
      vtkErrorMacro(<<"Could not extract pixeldata for " << (*dicom_files_p)[i].filename.c_str());
      }
    else // we have a dicom_pixeldata_pointer
      {
      if (temp_PixelRepresentation == 0) // "normal" unsigned pixeldata
        {
        if ((*si_iterator).BitsAllocated == 8)
          {
          for (int pidx = 0; pidx < numpixels_perslice; pidx++)
            {
            // HU = slope * SV + intercept
            // where: SV = actual unsigned int >> (highbit + 1 - bitsstored) && all valid stored bits on;
            // our output will be SIGNED shorts in this case...
            *(spixels + numpixels_perslice * i + pidx) =
              (short)INTROUND(temp_RescaleSlope * (double)((*(dicom_pixeldata_pointer_uc + pidx) >> rshift_factor) & validbits_mask) + temp_RescaleIntercept);
            }
          }
        else if (output->GetScalarType() == VTK_UNSIGNED_SHORT)// 16 BitsAllocated
          {
          // this means BitsStored is 16, so we output into unsigned short
          int pidx;
          for (pidx = 0; pidx < numpixels_perslice; pidx++)
            {
            // HU = slope * SV + intercept
            // where: SV = actual unsigned int >> (highbit + 1 - bitsstored) && all valid stored bits on;
            *(uspixels + numpixels_perslice * i + pidx) =
              (unsigned short)INTROUND(temp_RescaleSlope * (double)((*(dicom_pixeldata_pointer_us + pidx) >> rshift_factor) & validbits_mask) + temp_RescaleIntercept);
            }
          }
        else
          {
          // this means BitsStored is LESS than 16, so even although the pixeldata is unsigned,
          // we output into signed shorts, as the Rescale thingies can turn us into signed
          int pidx;
          for (pidx = 0; pidx < numpixels_perslice; pidx++)
            {
            // HU = slope * SV + intercept
            // where: SV = actual unsigned int >> (highbit + 1 - bitsstored) && all valid stored bits on;
            *(spixels + numpixels_perslice * i + pidx) =
              (short)INTROUND(temp_RescaleSlope * (double)((*(dicom_pixeldata_pointer_us + pidx) >> rshift_factor) & validbits_mask) + temp_RescaleIntercept);
            }
          }
        }
      else // two's complement
        {
        if ((*si_iterator).BitsAllocated == 8)
          {
          unsigned char temp_bits; 
          char* temp_bits_ptr;
          for (int pidx = 0; pidx < numpixels_perslice; pidx++)
            {
            temp_bits = (*(dicom_pixeldata_pointer_uc + pidx) >> rshift_factor) & validbits_mask;
            temp_bits_ptr = (char*)(&temp_bits); // here we get two's complement for free
            // and now we can get the HU again
            *(spixels + numpixels_perslice * i + pidx) = 
              (short)INTROUND(temp_RescaleSlope * (double)(*temp_bits_ptr) + temp_RescaleIntercept);
            }
          }
        else // 16 BitsAllocated
          {

          unsigned short temp_bits; 
          short* temp_bits_ptr;
          int pidx;
          for (pidx = 0; pidx < numpixels_perslice; pidx++)
            {
            temp_bits = (*(dicom_pixeldata_pointer_us + pidx) >> rshift_factor) & validbits_mask;
            temp_bits_ptr = (short*)(&temp_bits); // here we get two's complement for free
            // and now we can get the HU again
            *(spixels + numpixels_perslice * i + pidx) = 
              (short)INTROUND(temp_RescaleSlope * (double)(*temp_bits_ptr) + temp_RescaleIntercept);
            } // for (int pidx = 0; ...
          } // else 16 BitsAllocated
        } // else two's complement
      } // else we have a dicom_pixeldata_pointer

    // we really won't be using this again, so we delete all
    // fileformats so that we can have our memory back
    // (for years and years, this reader was using double the memory
    // it actually required)
    if ((*dicom_files_p)[i].fileformat)
      {
      delete (*dicom_files_p)[i].fileformat;
      (*dicom_files_p)[i].fileformat = NULL;
      }
    
    } // for (unsigned i = 0; i < dicom_files_p->size(); i++) ...
}

// ##########################################################################
unsigned long vtkDICOMVolumeReader::GetMTime(void)
{
  // this is what our parent would want
  unsigned long mtime = this->vtkStructuredPointsSource::GetMTime();

  // compare dicom_filenames_buffer and dicom_filenames
  bool vectors_equal = (dicom_filenames == dicom_filenames_buffer);

  if (!vectors_equal)
    {
    // the vectors differ, which means we are modified
    // 1. copy the vectors
    dicom_filenames = dicom_filenames_buffer;
    // 2. set the new modified time (there must be a less round-about way to do this!)
    vtkTimeStamp* tts = vtkTimeStamp::New();
    tts->Modified();
    mtime = tts->GetMTime();
    }

  return mtime;
}

std::list<series_instance>::iterator vtkDICOMVolumeReader::find_si_iterator(int idx)
{
  // the following depends on the currently selected SeriesInstanceIdx
  int si_cd = 0;
  std::list<series_instance>::iterator si_iterator;

  for (si_iterator = series_instances.begin(); si_cd != SeriesInstanceIdx && si_iterator != series_instances.end(); si_iterator++)
    {
    si_cd++;
    }
  return si_iterator;
  // if the index wasn't found, si_iterator will be == series_instances.end()
}

// ##########################################################################
void vtkDICOMVolumeReader::add_dicom_filename(const char* dicom_filename)
{
  // only the buffer changes, so we do not update MTime yet... only when
  // GetMTime() is called, we check for differences between this list
  // and dicom_filenames.  If there's no difference, MTime is NOT adapted
  dicom_filenames_buffer.push_back(std::string(dicom_filename));
}

// ##########################################################################
void vtkDICOMVolumeReader::clear_dicom_filenames(void)
{
  // we only clear the buffer list, not the internal one.  
  this->dicom_filenames_buffer.clear();
}


// ##########################################################################
int vtkDICOMVolumeReader::get_number_of_dicom_filenames(void)
{
  return dicom_filenames_buffer.size();
}

// ##########################################################################
const char *vtkDICOMVolumeReader::get_dicom_filename(int idx)
{
  if (idx < dicom_filenames_buffer.size())
    {
    return dicom_filenames_buffer[idx].c_str();
    }
  else
    {
    return NULL;
    }
}

// ##########################################################################
const char* vtkDICOMVolumeReader::GetSeriesInstanceUID(void)
{
  // we have to make sure that our data-structures are up-to-date
  this->UpdateInformation();
  std::list<series_instance>::iterator si_iterator = this->find_si_iterator(
    this->SeriesInstanceIdx);
  
  if (si_iterator == this->series_instances.end())
    {
    return NULL;
    }
  else
    {
    return (*si_iterator).SeriesInstanceUID.c_str();
    }
}

// ##########################################################################
const char *vtkDICOMVolumeReader::GetStudyDescription(void)
{
  this->UpdateInformation();
  std::list<series_instance>::iterator si_iterator = this->find_si_iterator(
    this->SeriesInstanceIdx);
  
  if (si_iterator == this->series_instances.end())
    {
    return NULL;
    }
  else
    {
    return (*si_iterator).misc_metadata.StudyDescription.c_str();
    }
}

// ##########################################################################
const char *vtkDICOMVolumeReader::GetReferringPhysician(void)
{
  this->UpdateInformation();
  std::list<series_instance>::iterator si_iterator = this->find_si_iterator(
    this->SeriesInstanceIdx);
  
  if (si_iterator == this->series_instances.end())
    {
    return NULL;
    }
  else
    {
    return (*si_iterator).misc_metadata.ReferringPhysician.c_str();
    }
}

// ##########################################################################
int vtkDICOMVolumeReader::GetMaximumSeriesInstanceIdx(void)
{
  this->UpdateInformation();
  return series_instances.size() - 1;
}


static char const rcsid[] =
"$Id$";

const char *vtkDICOMVolumeReader_rcsid(void)
{
  return rcsid;
}
