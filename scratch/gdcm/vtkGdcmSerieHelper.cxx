#include "vtkGdcmSerieHelper.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkGdcmSerieHelper, "$Revision: 1.38 $");

vtkStandardNewMacro(vtkGdcmSerieHelper);

vtkGdcmSerieHelper::vtkGdcmSerieHelper()
{
  this->helper = new gdcm::SerieHelper();
}

vtkGdcmSerieHelper::~vtkGdcmSerieHelper()
{
  delete this->helper;
}

void vtkGdcmSerieHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkGdcmSerieHelper::SetDirectory(char *dir_name)
{
  this->helper->SetDirectory(dir_name);
}

gdcm::FileList *vtkGdcmSerieHelper::GetFileList(char *serieUID)
{
  gdcm::FileList *cur_file_set = this->helper->GetFirstSingleSerieUIDFileSet();
  std::string cur_uid = this->helper->GetCurrentSerieUIDFileSetUID();

  while (cur_uid != serieUID && cur_file_set != NULL)
    {
    cur_file_set = this->helper->GetNextSingleSerieUIDFileSet();
    cur_uid = this->helper->GetCurrentSerieUIDFileSetUID();
    }

  if (cur_uid == serieUID)
    return cur_file_set;

  else
    return NULL;
}

int vtkGdcmSerieHelper::GetNumberOfUIDs(void)
{
  int total = 0;

  gdcm::FileList *cur_file_set = this->helper->GetFirstSingleSerieUIDFileSet();
  if (cur_file_set)
    total += 1;

  while (cur_file_set != NULL)
    {
    cur_file_set = this->helper->GetNextSingleSerieUIDFileSet();
    if (cur_file_set)
      total += 1;
    }

  return total;
}

char *vtkGdcmSerieHelper::GetUID(int idx)
{

  gdcm::FileList *cur_file_set = this->helper->GetFirstSingleSerieUIDFileSet();
  std::string cur_uid = this->helper->GetCurrentSerieUIDFileSetUID();
  int cur_idx = 0;

  while (cur_idx != idx && cur_file_set != NULL)
    {
    cur_file_set = this->helper->GetNextSingleSerieUIDFileSet();
    cur_idx += 1;
    }

  if (cur_idx == idx && cur_file_set)
    {
    strncpy(this->temp_uid, cur_uid.c_str(), 8191);
    return this->temp_uid;
    }
  else
    return NULL;
}

// void vtkGdcmSerieHelper::SetFileListOnReader(vtkGdcmReader *reader,
// 					     char *serieUID)
// {
//   gdcm::FileList *file_list = this->GetFileList(serieUID);
//   if (file_list)
//     {
//     reader->SetCoherentFileList(file_list);
//     }
//   else
//     {
//     vtkErrorMacro(<<"No filenames found for serie " << serieUID << ".");
//     }
// }
