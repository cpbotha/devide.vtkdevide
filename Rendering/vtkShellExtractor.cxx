// vtkShellExtractor.h copyright (c) 2002 by Charl P. Botha http://cpbotha.net/
// $Id: vtkShellExtractor.cxx,v 1.3 2003/05/06 11:34:47 cpbotha Exp $
// vtk class for extracting Udupa Shells

#include <math.h>
#include "vtkMatrix4x4.h"
#include "vtkObject.h"
#include "vtkObjectFactory.h"
#include "vtkShellExtractor.h"
#include "vtkTransform.h"

#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkLongArray.h"
#include "vtkPointData.h"
#include "vtkUnsignedLongArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"

vtkStandardNewMacro(vtkShellExtractor);

#include <vector>
using namespace std;

/**
 * Template class for doing the actual work behind shell extraction.
 * 
 * @param data_ptr Pointer to xyz volume that the shell must be extracted from.
 * @param xlen Number of voxels in x-direction.
 * @param ylen Number of voxels in y-direction.
 * @param zlen Number of voxels in z-direction.
 * @param x_spacing Object space distance between voxels in X.
 * @param y_spacing Object space distance between voxels in Y.
 * @param z_spacing Object space distance between voxels in Y.
 * @param OpacityTF The opacity function that you're rendering with.
 * @param OmegaL Lower opacity limit.
 * @param OmegaH Upper opacity limit.
 * @param vectorD Pointer to vector of ShellVoxels which will be filled with
 * the result of the extraction.
 * @param P array of ints of size ylen*zlen
 */
template <class T> 
static void ExtractShell(T* data_ptr,
                         vtkImageData* Input,
                         vtkImageData* GradientImageData,
                         vtkPiecewiseFunction* OpacityTF,
                         vtkColorTransferFunction* ColourTF,
                         float OmegaL, float OmegaH, float min_gradmag,
                         vector<ShellVoxel>* vectorD, int P[])
{
    T* dptr = data_ptr;
    T nbrs[6]; // the six neighbours
    float gnbrs[6]; // the six neighbours in the gradient volume
    float nbr_os[6]; // opacities of the six neighbours
    float t_gradmag; // temp gradient magnitude

    int nbrv_lt_oh; // flag variable: neighbourhood voxel < OmegaH ?
    float temp_rgb[3];

    // important variables
    int xlen = Input->GetDimensions()[0];
    int ylen = Input->GetDimensions()[1];
    int zlen = Input->GetDimensions()[2];
    float x_spacing = Input->GetSpacing()[0];
    float y_spacing = Input->GetSpacing()[1];
    float z_spacing = Input->GetSpacing()[2];
    float x_orig = Input->GetOrigin()[0];
    float y_orig = Input->GetOrigin()[1];
    float z_orig = Input->GetOrigin()[2];   

    //vector<ShellVoxel> vectorD;
    // initialise the D vector
    vectorD->clear();
    ShellVoxel temp_sv;

    // initialise P to all "no shell voxels in these rows (X)"
    memset(P, -1, ylen * zlen * sizeof(int) * 2);
    int Pidx, prevPidx;

    int xstep = 1;
    int ystep = xlen;
    int zstep = xlen * ylen;

    float dxm2 = 2.0 * x_spacing;
    float dym2 = 2.0 * y_spacing;
    float dzm2 = 2.0 * z_spacing;

    vtkTransform *voxelsToVolumeTransform = vtkTransform::New();
    voxelsToVolumeTransform->Identity();
    voxelsToVolumeTransform->Translate(Input->GetOrigin());
    voxelsToVolumeTransform->Scale(Input->GetSpacing());
    vtkMatrix4x4 *voxelsToVolumeMatrix = vtkMatrix4x4::New();
    voxelsToVolumeMatrix->DeepCopy(voxelsToVolumeTransform->GetMatrix());
    voxelsToVolumeTransform->Delete();

    for (int z = 0; z < zlen; z++)
    {
        for (int y = 0; y < ylen; y++)
        {
            for (int x = 0; x < xlen; x++)
            {
                temp_sv.Value = (float)(*dptr);
                // look the suxor up
                temp_sv.Opacity = OpacityTF->GetValue(temp_sv.Value);
                // first check if it's opaque enough
                if (temp_sv.Opacity > OmegaL)
                {
                    // initially, we assume none of the neighbouring voxels
                    // is < OmegaH
                    nbrv_lt_oh = 0;
                    // and zero the neighbour code
                    temp_sv.nbrOCode = 0;

                    // cache the six neighbour values since we'll probably
                    // be needing them again; the cpu cache doesn't like the
                    // read we're doing here... :(
                    // note that we're zero-padding for both the normal
                    // calculation and neighbourhood thingies

                    // X ************************************************
                    if (x > 0)
                    {
                        nbrs[0] = *(dptr - xstep);
                        nbr_os[0] = OpacityTF->GetValue((float)(nbrs[0]));
                        if (nbr_os[0] >= OmegaH)
                            temp_sv.nbrOCode |= 0x1; // activate bit 0
                        if (nbr_os[0] <= OmegaH)
                            nbrv_lt_oh = 1;
                    }
                    else
                    {
                        nbrs[0] = 0;
                        nbr_os[0] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }

                    if (x < xlen - 1)
                    {
                        nbrs[1] = *(dptr + xstep);
                        nbr_os[1] = OpacityTF->GetValue((float)(nbrs[1]));
                        if (nbr_os[1] >= OmegaH)
                            temp_sv.nbrOCode |= 0x2; // activate bit 1
                        if (nbr_os[1] <= OmegaH)
                            nbrv_lt_oh = 1;
                    }
                    else
                    {
                        nbrs[1] = 0;
                        nbr_os[1] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }



                    // Y ************************************************
                    if (y > 0)
                    {
                        nbrs[2] = *(dptr - ystep);
                        nbr_os[2] = OpacityTF->GetValue((float)(nbrs[2]));
                        if (nbr_os[2] >= OmegaH)
                            temp_sv.nbrOCode |= 0x4; // activate bit 2
                        if (nbr_os[2] <= OmegaH)
                            nbrv_lt_oh = 1;
                    }
                    else
                    {
                        nbrs[2] = 0;
                        nbr_os[2] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }


                    if (y < ylen - 1)
                    {
                        nbrs[3] = *(dptr + ystep);
                        nbr_os[3] = OpacityTF->GetValue((float)(nbrs[3]));
                        if (nbr_os[3] >= OmegaH)
                            temp_sv.nbrOCode |= 0x8; // activate bit 3
                        if (nbr_os[3] <= OmegaH)
                            nbrv_lt_oh = 1;
                    }
                    else
                    {
                        nbrs[3] = 0;
                        nbr_os[3] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }


                    // Z ************************************************		   
                    if (z > 0)
                    {
                        nbrs[4] = *(dptr - zstep);
                        nbr_os[4] = OpacityTF->GetValue((float)(nbrs[4]));
                        if (nbr_os[4] >= OmegaH)
                            temp_sv.nbrOCode |= 0x10; // activate bit 4
                        if (nbr_os[4] <= OmegaH)
                            nbrv_lt_oh = 1;

                    }
                    else
                    {
                        nbrs[4] = 0;
                        nbr_os[4] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }


                    if (z < zlen - 1)
                    {
                        nbrs[5] = *(dptr + zstep);
                        nbr_os[5] = OpacityTF->GetValue((float)(nbrs[5]));
                        if (nbr_os[5] >= OmegaH)
                            temp_sv.nbrOCode |= 0x20; // activate bit 5
                        if (nbr_os[5] <= OmegaH)
                            nbrv_lt_oh = 1;
                    }
                    else
                    {
                        nbrs[5] = 0;
                        nbr_os[5] = 0.0;
                        // also, since 0.0 < OmegaL <= OmegaH < 1.0, we have a
                        // neighbour with opacity <= OmegaH
                        nbrv_lt_oh = 1;
                        // we don't have to adapt the neighbour opacity code,
                        // since it's zeroed by default
                    }


                    // DONE with neighbour caching, lookups and neighbour code

                    // only if one of the neighbouring voxels was <= OmegaH
                    // is this a valid shell voxel
                    if (nbrv_lt_oh)
                        //if (1)
                    {
                        if (!GradientImageData)
                        {
                            // calculate normal
                            // NOTE: we are dividing by 2.0 * spacing as opengl wants
                            // these things to be NORMALISED in world space (although
                            // we'll be rendering in voxel space)
                            temp_sv.Normal[0] = ((float)nbrs[0] - (float)nbrs[1]) / dxm2;
                            temp_sv.Normal[1] = ((float)nbrs[2] - (float)nbrs[3]) / dxm2;
                            temp_sv.Normal[2] = ((float)nbrs[4] - (float)nbrs[5]) / dxm2;
                        }
                        else
                        {
                            // make use of the GradientImageVolume for calculating gradients
                            if (x > 0)
                                gnbrs[0] = GradientImageData->GetScalarComponentAsFloat(x - 1, y, z, 0);
                            else
                                gnbrs[0] = 0.0;
                            
                            if (x < xlen - 1)
                                gnbrs[1] = GradientImageData->GetScalarComponentAsFloat(x + 1, y, z, 0);
                            else
                                gnbrs[1] = 0.0;

                            if (y > 0)
                                gnbrs[2] = GradientImageData->GetScalarComponentAsFloat(x, y - 1, z, 0);
                            else
                                gnbrs[2] = 0.0;
                            
                            if (y < ylen - 1)
                                gnbrs[3] = GradientImageData->GetScalarComponentAsFloat(x, y + 1, z, 0);
                            else
                                gnbrs[3] = 0.0;

                            if (z > 0)
                                gnbrs[4] = GradientImageData->GetScalarComponentAsFloat(x, y, z - 1, 0);
                            else
                                gnbrs[4] = 0.0;

                            if (z < zlen - 1)
                                gnbrs[5] = GradientImageData->GetScalarComponentAsFloat(x, y, z + 1, 0);
                            else
                                gnbrs[5] = 0.0;

                            temp_sv.Normal[0] = ((float)gnbrs[0] - (float)gnbrs[1]) / dxm2;
                            temp_sv.Normal[1] = ((float)gnbrs[2] - (float)gnbrs[3]) / dxm2;
                            temp_sv.Normal[2] = ((float)gnbrs[4] - (float)gnbrs[5]) / dxm2;
                        }

                        t_gradmag = sqrt(temp_sv.Normal[0] * temp_sv.Normal[0] +
                                         temp_sv.Normal[1] * temp_sv.Normal[1] +
                                         temp_sv.Normal[2] * temp_sv.Normal[2]);

                        if (t_gradmag > min_gradmag)
                        {
                            temp_sv.Normal[0] /= t_gradmag;
                            temp_sv.Normal[1] /= t_gradmag;
                            temp_sv.Normal[2] /= t_gradmag;
                        }
                        else
                        {
                            // the normal is too small, shouldn't count
                            temp_sv.Normal[0] = temp_sv.Normal[1] = temp_sv.Normal[2] = 0.0;
                        }


                        // now we can fill up that structure
                        // Value is done
                        // Opacity is done
                        // Normal is done
                        // nbrOCode is done
                        temp_sv.x = x; // set integer x
                        // and update the voxel coords in voxel space
                        temp_sv.volCoords[0] = x_orig + (float)x * x_spacing;
                        temp_sv.volCoords[1] = y_orig + (float)y * y_spacing;
                        temp_sv.volCoords[2] = z_orig + (float)z * z_spacing;                       

                        if (!ColourTF)
                            temp_sv.Red = temp_sv.Green = temp_sv.Blue = 1.0;
                        else
                        {
                            ColourTF->GetColor(temp_sv.Value, temp_rgb);
                            temp_sv.Red = temp_rgb[0];
                            temp_sv.Green = temp_rgb[1];
                            temp_sv.Blue = temp_rgb[2];
                        }

                        vectorD->push_back(temp_sv);
                        Pidx = (z*ylen + y) * 2;
                        if (P[Pidx] == -1)
                        {
                            // this is the first voxel in this x-row, set its
                            // index in the P structure
                            P[Pidx] = vectorD->size() - 1;
                            // this means (by definition) that the previous row
                            // is done (or had no shell voxels at all, in which
                            // case it has no length, and in which case we
                            // have to search farther back for a row to tally
                            // up)

                            // we have to travel back to find the previous
                            // row with shell voxels in it
                            prevPidx = Pidx - 2;
                            while (prevPidx >= 0 && P[prevPidx] == -1)
                                prevPidx-=2;
                            // if we found a valid prevPidx, we can tally,
                            // if not, Pidx is the FIRST row
                            if (prevPidx >= 0)
                            {
                                P[prevPidx+1] = P[Pidx] - P[prevPidx];
                            }
                        }
                    } // if (nbrv_lt_oh) ...
                } // if (temp_sv.Opacity > OmegaL ...

                // make sure to increment dptr!!!
                dptr++;

            } // for (int x = 0 ...
        } // for (int y = 0 ...
        cout << "row " << z << " done." << endl;
    } // for (int z = 0 ...

    // we're done, yahooooo!
    // do some cleanup
    voxelsToVolumeMatrix->Delete();
}


vtkShellExtractor::vtkShellExtractor()
{
    this->Input = NULL;
    this->OmegaL = 0.0;
    this->OmegaH = 1.0;
    this->OpacityTF = NULL;
    this->ColourTF = NULL;
    this->GradientImageData = NULL;
    this->D = NULL;
    this->P = NULL;
}

vtkShellExtractor::~vtkShellExtractor()
{
    this->SetInput(NULL);
    this->SetOpacityTF(NULL);
    this->SetColourTF(NULL);
    this->SetGradientImageData(NULL);
    if (this->D)
        delete this->D;
    if (this->P)
        delete this->P;
}

void vtkShellExtractor::Update(void)
{
    if (!this->Input)
    {
        vtkErrorMacro(<< "No input set to calculate shell on.");
        return;
    }

    if (!this->OpacityTF)
    {
        vtkErrorMacro(<< "No opacity transfer function set.");
        return;
    }

    if (this->GetMTime() > this->ExtractTime ||
        this->OpacityTF->GetMTime() > this->ExtractTime ||
        this->ColourTF && this->ColourTF->GetMTime() > this->ExtractTime ||
        this->Input->GetMTime() > this->ExtractTime ||
        this->GradientImageData && 
        this->GradientImageData->GetMTime() > this->ExtractTime ||
        !this->D)
    {
        // make sure our input is up to date
        this->Input->UpdateInformation();
        this->Input->SetUpdateExtentToWholeExtent();
        this->Input->Update();

        // and the transfer function too
        this->OpacityTF->Update();


        // we can only use the GradientImageData if it has exactly the same number
        // of elements as the Input
        vtkImageData* TempGradientImageData = NULL;

        // update the GradientImageData IF we have it
        if (this->GradientImageData)
        {
            this->GradientImageData->UpdateInformation();
            this->GradientImageData->SetUpdateExtentToWholeExtent();
            this->GradientImageData->Update();

            if (this->GradientImageData->GetDimensions()[0] == this->Input->GetDimensions()[0] &&
                this->GradientImageData->GetDimensions()[1] == this->Input->GetDimensions()[1] &&
                this->GradientImageData->GetDimensions()[2] == this->Input->GetDimensions()[2])
            {
                TempGradientImageData = GradientImageData;
            }
            else
            {
                vtkWarningMacro("<< GradientImageData and Input have different dimensions.");
            }
            
        }

        // we're going to iterate through each voxel in the volume and:
        // 1. if the voxel opacity is > OmegaL and it has at least one voxel 
        //    in its neighbourhood that is < OmegaH then
        //    a. calculate object space normal and store
        //    b. extract neighbourhood code
        //    c. store everything and append to D and update P

        int idims[3];
        this->Input->GetDimensions(idims);
        float ispacing[3];
        this->Input->GetSpacing(ispacing);

        // do stuff here
        vtkDataArray* scalars = this->Input->GetPointData()->GetScalars();
        if (scalars == NULL)
        {
            vtkErrorMacro(<< "No scalars to extract shell from!");
            return;
        }

        if (this->P)
            delete this->P;
        // need array of Y * Z ints
        this->P = new int[idims[1] * idims[2] * 2];

        vector<ShellVoxel> vectorD;
        float min_gradmag = 0.0;

        switch (scalars->GetDataType())
        {
        case VTK_CHAR:
        {
            char* data_ptr = ((vtkCharArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_UNSIGNED_CHAR:
        {
            unsigned char* data_ptr = ((vtkUnsignedCharArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_SHORT:
        {
            short* data_ptr = ((vtkShortArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_UNSIGNED_SHORT:
        {
            unsigned short* data_ptr = ((vtkUnsignedShortArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_INT:
        {
            int* data_ptr = ((vtkIntArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_UNSIGNED_INT:
        {
            unsigned int* data_ptr = ((vtkUnsignedIntArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_LONG:
        {
            long* data_ptr = ((vtkLongArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_UNSIGNED_LONG:
        {
            unsigned long* data_ptr = ((vtkUnsignedLongArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_FLOAT:
        {
            float* data_ptr = ((vtkFloatArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        case VTK_DOUBLE:
        {
            double* data_ptr = ((vtkDoubleArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorD, this->P);
        }
        break;
        default:
            vtkErrorMacro(<< "Unknown type of scalar data!");
        }

        // now we should copy vectorD over to D
        if (this->D)
            delete this->D;

        this->D = new ShellVoxel[vectorD.size()];
        for (int i = 0; i < vectorD.size(); i++)
            memcpy(this->D + i, &(vectorD[i]), sizeof(ShellVoxel));

        cout << vectorD.size() << " shell voxels found." << endl;

        // stop doing stuff
        this->ExtractTime.Modified();
    }
}

const int vtkShellExtractor::shell_noc_visibility_lut[8][64] =
{
    // octant 0, bits 0 & 2 & 4
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0
    },
    // octant 1, bits 1 & 2 & 4
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0
    },
    // octant 2, bits 0 & 3 & 4
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0
    },
    // octant 3, bits 1 & 3 & 4
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0
    },
    // octant 4, bits 0 & 2 & 5
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0,
     1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0
    },
    // octant 5, bits 1 & 2 & 5
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
     1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0
    },
    // octant 6, bits 0 & 3 & 5
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0
    },
    // octant 7, bits 1 & 3 & 5
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0
    }
};

