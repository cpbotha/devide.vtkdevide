// vtkHDFVolumeReader.cxx copyright 2001,2002 Charl P. Botha <cpbotha@ieee.org>
// $Id: vtkHDFVolumeReader.cxx,v 1.2 2004/01/15 11:00:55 cpbotha Exp $
// vtk class for reading scalar volume data

#include "vtkHDFVolumeReader.h"
#include <vtkObjectFactory.h>
#include <vtkStructuredPoints.h>
#include <vtkUnsignedShortArray.h>
#include <string.h>

#define DATA_SPACING_ATTR_NAME "data_spacing"

vtkHDFVolumeReader::vtkHDFVolumeReader()
{
   // we have no valid spacing yet, so we show it this way
   DataSpacing[0] = DataSpacing[1] = DataSpacing[2] = 0.0;
   DataOrigin[0] = DataOrigin[1] = DataOrigin[2] = 0.0;
   data_dimensions[0] = data_dimensions[1] = data_dimensions[2] = 0;
   sd_id = sds_id = -1;
   FileName = NULL;
}

vtkHDFVolumeReader::~vtkHDFVolumeReader()
{
   close_hdf();
   // deallocate space allocated for this by SetStringMacro
   if (FileName)
   {
      delete[] FileName;
   }
}

int vtkHDFVolumeReader::open_hdf(void)
{
   // make sure we don't have any open files lying around (and this
   // well set status_OK to true)
   close_hdf();

   int32 status;

   if (!FileName || strlen(FileName)==0)
   {
      vtkErrorMacro(<<"FileName not set");
      return 0;
   }

   if ((sd_id = SDstart(FileName, DFACC_READ)) == FAIL)
   {
      vtkErrorMacro(<<"Could not open HDF file.");
      return 0;
   }

   if ((sds_id = SDselect(sd_id, 0)) == FAIL)
   {
      vtkErrorMacro(<<"Could not get selector for SDS in open HDF.");
      return 0;
   }

   return 1;
}

void vtkHDFVolumeReader::close_hdf(void)
{
   if (sd_id >= 0)
   {
      // if we haven't de-initialised yet, do it now
      SDendaccess(sds_id);
      SDend(sd_id);
      // and other methods should know that we're closed
      sd_id = -1;
   }
}

void vtkHDFVolumeReader::ExecuteInformation(void)
{
    if (!open_hdf())
    {
        return;
   }

   char name[MAX_NC_NAME];
   int32 dim_sizes[MAX_VAR_DIMS];
   int32 rank, data_type, n_attrs; // number of dimensions, we want 3
   int32 status;

   // FIXME: check here that the data is 3 dimensional and consists of 
   // unsigned ints
   status = SDgetinfo(sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);

   if (status == FAIL)
   {
      vtkErrorMacro(<<"Error retrieving SDS metadata.");
      return;
   }

   // now get our devide1 specific attributes (we should abstract this and  move it out)
   float data_spacing[3];
   bool attr_error = false;
   int attr_idx = SDfindattr(sds_id, DATA_SPACING_ATTR_NAME);
   if (attr_idx >= 0)
   {
      int32 data_type, n_values;
      char attr_name[MAX_NC_NAME];
      status = SDattrinfo(sds_id, attr_idx, attr_name, &data_type, &n_values);
      if (status != FAIL && n_values == 3 && data_type == DFNT_FLOAT32)
         status = SDreadattr(sds_id, attr_idx, data_spacing);
      else
         attr_error = true;
   }
   else
      attr_error = true;
   // if we were successful in reading the spacing attribute, we can set this up
   if (!attr_error)
   {
      DataSpacing[0] = data_spacing[0];
      DataSpacing[1] = data_spacing[1];
      DataSpacing[2] = data_spacing[2];
   }
   else
   {
      DataSpacing[0] = DataSpacing[1] = DataSpacing[2] = 1.0;
   }

   // HDF has the slowest changing dimension first, we want it the other way round
   data_dimensions[0] = dim_sizes[2];
   data_dimensions[1] = dim_sizes[1];
   data_dimensions[2] = dim_sizes[0];

   // Transfer all information to output ---------------------------------
   // --------------------------------------------------------------------------
   
   // get pointer to our output
   vtkStructuredPoints* output = this->GetOutput();

   output->SetDimensions(data_dimensions);
   output->SetSpacing(DataSpacing);
   output->SetOrigin(DataOrigin);
   output->SetWholeExtent(0, data_dimensions[0]-1, 0, data_dimensions[1]-1, 0, data_dimensions[2]-1);

   output->SetScalarType(VTK_UNSIGNED_SHORT);
   output->SetNumberOfScalarComponents(1);
}

void vtkHDFVolumeReader::ExecuteData(vtkDataObject* out)
{
   if (!open_hdf())
   {
      return;
   }
   
   if (sd_id < 0)
   {
      vtkErrorMacro(<<"vtkHDFVolumeReader::ExecuteData() called un-initialised");
      return;
   }

   vtkStructuredPoints* output = vtkStructuredPoints::SafeDownCast(out);
   if (!output)
   {
      vtkWarningMacro("ExecuteData with non-vtkStructuredPoints output!");
      return;
   }
  
   // we're going to read the whole dataset, so make sure the extent is good
   // we could just do SetExtent(output->GetUpdateExtent()), but then we would have to
   // be more clever with our actual reading routines, i.e. only read the data necessary
   // to create the UpdateExtent in the output
   output->SetExtent(output->GetWholeExtent());
   // AllocateScalars will make use of the information set on the output by 
   // ExecuteInformation() to perform the necessary allocation
   output->AllocateScalars();

   // find a pointer that we can use
   unsigned short* pixels = (unsigned short*)output->GetScalarPointer();

   // fire up the HDF and read the data
   int32 old_dims[3];
   old_dims[0] = data_dimensions[2];
   old_dims[1] = data_dimensions[1];
   old_dims[2] = data_dimensions[0];
   int32 start[3];
   start[0] = start[1] = start[2] = 0;

   int status = SDreaddata(sds_id, start, NULL, old_dims, pixels);
   if (status == FAIL)
   {
      vtkErrorMacro(<<"Unable to read dataset.");
      return;
   }
}

vtkHDFVolumeReader* vtkHDFVolumeReader::New()
{
   vtkObject* ret = vtkObjectFactory::CreateInstance("vtkHDFVolumeReader");
   if (ret)
   {
      return(vtkHDFVolumeReader*)ret;
   }
   return new vtkHDFVolumeReader;
}

static char const rcsid[] =
"$Id: vtkHDFVolumeReader.cxx,v 1.2 2004/01/15 11:00:55 cpbotha Exp $";

const char *vtkHDFVolumeReader_rcsid(void)
{
   return rcsid;
}
