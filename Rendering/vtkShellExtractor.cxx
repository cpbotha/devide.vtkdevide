// vtkShellExtractor.h copyright (c) 2003 
// by Charl P. Botha cpbotha@ieee.org 
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id: vtkShellExtractor.cxx,v 1.13 2004/06/29 22:12:51 cpbotha Exp $
// vtk class for extracting Udupa Shells

/*
 * This software is licensed exclusively for research use by Jorit Schaap.
 * Any modifications made to this software shall be sent to the author for 
 * possible inclusion in future versions.  Ownership and copyright of said 
 * modifications shall be ceded to the author.
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
			 int GradientImageIsGradient,
                         vtkPiecewiseFunction* OpacityTF,
                         vtkColorTransferFunction* ColourTF,
                         float OmegaL, float OmegaH, float min_gradmag,
                         vector<ShellVoxel>* vectorDx, int Px[],
                         vector<ShellVoxel>* vectorDy, int Py[],
                         vector<ShellVoxel>* vectorDz, int Pz[])
{
    T* dptr = data_ptr;
    T nbrs[6]; // the six neighbours
    float gnbrs[6]; // the six neighbours in the gradient volume
    float nbr_os[6]; // opacities of the six neighbours
    float t_gradmag; // temp gradient magnitude

    int nbrv_lt_oh; // flag variable: neighbourhood voxel < OmegaH ?
    double temp_rgb[3];

    // important variables
    int xlen = Input->GetDimensions()[0];
    int ylen = Input->GetDimensions()[1];
    int zlen = Input->GetDimensions()[2];
    double x_spacing = Input->GetSpacing()[0];
    double y_spacing = Input->GetSpacing()[1];
    double z_spacing = Input->GetSpacing()[2];
    double x_orig = Input->GetOrigin()[0];
    double y_orig = Input->GetOrigin()[1];
    double z_orig = Input->GetOrigin()[2];   

    //vector<ShellVoxel> vectorD;
    // initialise the D vector
    vectorDx->clear();
    vectorDy->clear();
    vectorDz->clear();
      
    ShellVoxel temp_sv;

    // initialise P to all "no shell voxels in these rows (X)"
    memset(Px, -1, ylen * zlen * sizeof(int) * 2);
    memset(Py, -1, ylen * zlen * sizeof(int) * 2);
    memset(Pz, -1, ylen * zlen * sizeof(int) * 2);
    
    int Pidx, prevPidx;

    int xstep = 1;
    int ystep = xlen;
    int zstep = xlen * ylen;

    double dxm2 = 2.0 * x_spacing;
    double dym2 = 2.0 * y_spacing;
    double dzm2 = 2.0 * z_spacing;
 
   double tempValue;

    for (int z = 0; z < zlen; z++)
    {
        for (int y = 0; y < ylen; y++)
        {
            for (int x = 0; x < xlen; x++)
            {
                //temp_sv.Value = (float)(*dptr);
                tempValue = (double)(*dptr);
                // look the suxor up
                temp_sv.Opacity = OpacityTF->GetValue(tempValue);
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
                        nbr_os[0] = OpacityTF->GetValue((double)(nbrs[0]));
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
                        nbr_os[1] = OpacityTF->GetValue((double)(nbrs[1]));
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
                        nbr_os[2] = OpacityTF->GetValue((double)(nbrs[2]));
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
                        nbr_os[3] = OpacityTF->GetValue((double)(nbrs[3]));
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
                        nbr_os[4] = OpacityTF->GetValue((double)(nbrs[4]));
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
                        nbr_os[5] = OpacityTF->GetValue((double)(nbrs[5]));
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
                            temp_sv.Normal[0] = ((double)nbrs[0] - (double)nbrs[1]) / dxm2;
                            temp_sv.Normal[1] = ((double)nbrs[2] - (double)nbrs[3]) / dym2;
                            temp_sv.Normal[2] = ((double)nbrs[4] - (double)nbrs[5]) / dzm2;
                        }
                        else if (!GradientImageIsGradient)
                        {
                            // make use of the GradientImageVolume for calculating gradients
                            if (x > 0)
                                gnbrs[0] = GradientImageData->GetScalarComponentAsDouble(x - 1, y, z, 0);
                            else
                                gnbrs[0] = 0.0;
                            
                            if (x < xlen - 1)
                                gnbrs[1] = GradientImageData->GetScalarComponentAsDouble(x + 1, y, z, 0);
                            else
                                gnbrs[1] = 0.0;

                            if (y > 0)
                                gnbrs[2] = GradientImageData->GetScalarComponentAsDouble(x, y - 1, z, 0);
                            else
                                gnbrs[2] = 0.0;
                            
                            if (y < ylen - 1)
                                gnbrs[3] = GradientImageData->GetScalarComponentAsDouble(x, y + 1, z, 0);
                            else
                                gnbrs[3] = 0.0;

                            if (z > 0)
                                gnbrs[4] = GradientImageData->GetScalarComponentAsDouble(x, y, z - 1, 0);
                            else
                                gnbrs[4] = 0.0;

                            if (z < zlen - 1)
                                gnbrs[5] = GradientImageData->GetScalarComponentAsDouble(x, y, z + 1, 0);
                            else
                                gnbrs[5] = 0.0;

                            temp_sv.Normal[0] = ((double)gnbrs[0] - (double)gnbrs[1]) / dxm2;
                            temp_sv.Normal[1] = ((double)gnbrs[2] - (double)gnbrs[3]) / dym2;
                            temp_sv.Normal[2] = ((double)gnbrs[4] - (double)gnbrs[5]) / dzm2;
                        }
			else
			{
			  temp_sv.Normal[0] = GradientImageData->GetScalarComponentAsDouble(x,y,z,0) / dxm2;
			  temp_sv.Normal[1] = GradientImageData->GetScalarComponentAsDouble(x,y,z,1) / dym2;
			  temp_sv.Normal[2] = GradientImageData->GetScalarComponentAsDouble(x,y,z,2) / dzm2;
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
                        //temp_sv.volCoords[0] = x_orig + (float)x * x_spacing;
                        //temp_sv.volCoords[1] = y_orig + (float)y * y_spacing;
                        //temp_sv.volCoords[2] = z_orig + (float)z * z_spacing;                       

                        if (!ColourTF)
                            temp_sv.Red = temp_sv.Green = temp_sv.Blue = 1.0;
                        else
                        {
                            ColourTF->GetColor(tempValue, temp_rgb);
                            temp_sv.Red = temp_rgb[0];
                            temp_sv.Green = temp_rgb[1];
                            temp_sv.Blue = temp_rgb[2];
                        }

                        vectorDx->push_back(temp_sv);
                        Pidx = (z*ylen + y) * 2;
                        if (Px[Pidx] == -1)
                        {
                            // this is the first voxel in this x-row, set its
                            // index in the P structure
                            Px[Pidx] = vectorDx->size() - 1;
                            // this means (by definition) that the previous row
                            // is done (or had no shell voxels at all, in which
                            // case it has no length, and in which case we
                            // have to search farther back for a row to tally
                            // up)

                            // we have to travel back to find the previous
                            // row with shell voxels in it
//                             prevPidx = Pidx - 2;
//                             while (prevPidx >= 0 && P[prevPidx] == -1)
//                                 prevPidx-=2;
//                             // if we found a valid prevPidx, we can tally,
//                             // if not, Pidx is the FIRST row
//                             if (prevPidx >= 0)
//                             {
//                                 P[prevPidx+1] = P[Pidx] - P[prevPidx];
//                             }
                        } // if (P[Pidx] == -1) ...
                    } // if (nbrv_lt_oh) ...
                } // if (temp_sv.Opacity > OmegaL ...

                // make sure to increment dptr!!!
                dptr++;

            } // for (int x = 0 ...

            // now make sure that we COMPLETE the current row
            Pidx = (z*ylen + y) * 2;
            if (Px[Pidx] != -1 && Px[Pidx+1] == -1)
              {
              // this means that the beginning of the row is indicated
              // but not the run length
              Px[Pidx + 1] = vectorDx->size() - Px[Pidx];
              }
            
        } // for (int y = 0 ...
        // FIXME: add progress event here (or something)
        //cout << "row " << z << " done." << endl;
    } // for (int z = 0 ...

    // we're done, yahooooo!
    // well not quite... let's build up the y and z Ps and Ds.

    // we need a temporary pointer matrix so we can charge through Px
    // and Dx to generate the others...
    ShellVoxel **pointerWallZY = new ShellVoxel*[zdim * ydim];
    ShellVoxel* tempDptr;
    int tempDoffset;
    
    // P is Z * X
    // 1. initialize pointerWallZY for the x -> y run
    memset(pointerWallZY, 0, zdim * ydim * sizeof(ShellVoxel *));
    for (int z = 0; z < zlen; z++)
      {
      for (int y = 0; y < ylen; y++)
        {
        // for each position in P we have an offset and a length!
        tempDoffset = Px[(z*ydim + y) * 2];
        if (tempDoffset >= 0)
          {
          pointerWall[z*ydim + y] = vectorDx(tempDoffset);
          }
        }
      }
    
    for (int z = 0; z < zlen; z++)
      {
      for (int x = 0; x < xlen; x++)
        {
        for (int y = 0; y < ylen; y++)
          {

          // find the Dptr for this x,y,z if it exists
          // arghhh... need maybe P-wall with also run-lengths so that
          // we know when to stop for a particular z,y combo
          // hmmmm, copy of P. at each successful Dptr access, we
          // decrement the run-length.  When the run-length reaches
          // zero, we're done.  The offset is obviously increased
          // every time
          Pidx = (z * ydim + y) * 2;
          if (pWallZY[Pidx] >= 0)
            {
            // find and access, then decrement run-length
            }
          
          
          } // for (int y = 0 ...
        } // for (int x = 0 ...
      } // for (int z = 0 ...
}


vtkShellExtractor::vtkShellExtractor()
{
    this->Input = NULL;
    this->OmegaL = 0.0;
    this->OmegaH = 1.0;
    this->OpacityTF = NULL;
    this->ColourTF = NULL;
    this->GradientImageData = NULL;
    this->GradientImageIsGradient = 0;
    this->Dx = this->Dy = this->Dz = NULL;
    this->Px = this->Py = this->Pz = NULL;
}

vtkShellExtractor::~vtkShellExtractor()
{
    this->SetInput(NULL);
    this->SetOpacityTF(NULL);
    this->SetColourTF(NULL);
    this->SetGradientImageData(NULL);
    
    if (this->Dx)
      delete this->Dx;
    
    if (this->Dy)
      delete this->Dy;

    if (this->Dz)
      delete this->Dz;
    
    if (this->Px)
      delete this->Px;

    if (this->Py)
      delete this->Py;

    if (this->Pz)
      delete this->Pz;
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
        !this->Dx)
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

            if (this->GradientImageIsGradient && 
                this->GradientImageData->GetNumberOfScalarComponents() < 3)
            {
                vtkWarningMacro("<< GradientImageIsData ON, but < 3 components in GradientImage.  Setting to OFF.");
                this->GradientImageIsGradient = 0;
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
        double ispacing[3];
        this->Input->GetSpacing(ispacing);

        // do stuff here
        vtkDataArray* scalars = this->Input->GetPointData()->GetScalars();
        if (scalars == NULL)
        {
            vtkErrorMacro(<< "No scalars to extract shell from!");
            return;
        }

        if (this->Px)
          delete this->Px;

        if (this->Py)
          delete this->Py;

        if (this->Pz)
          delete this->Pz;
        
        // need array of Z * Y ints
        this->Px = new int[idims[2] * idims[1] * 2];
        // Z * X ints
        this->Py = new int[idims[2] * idims[0] * 2];
        // Y * X ints
        this->Pz = new int[idims[1] * idims[0] * 2];

        vector<ShellVoxel> vectorDx, vectorDy, vectorDz;
        double min_gradmag = 0.0;

        switch (scalars->GetDataType())
        {
        case VTK_CHAR:
        {
            char* data_ptr = ((vtkCharArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_UNSIGNED_CHAR:
        {
            unsigned char* data_ptr = ((vtkUnsignedCharArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_SHORT:
        {
            short* data_ptr = ((vtkShortArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);

        }
        break;
        case VTK_UNSIGNED_SHORT:
        {
            unsigned short* data_ptr = ((vtkUnsignedShortArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_INT:
        {
            int* data_ptr = ((vtkIntArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_UNSIGNED_INT:
        {
            unsigned int* data_ptr = ((vtkUnsignedIntArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData, 
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_LONG:
        {
            long* data_ptr = ((vtkLongArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_UNSIGNED_LONG:
        {
            unsigned long* data_ptr = ((vtkUnsignedLongArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_FLOAT:
        {
            float* data_ptr = ((vtkFloatArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        case VTK_DOUBLE:
        {
            double* data_ptr = ((vtkDoubleArray*)scalars)->GetPointer(0);
            ExtractShell(data_ptr, this->Input, TempGradientImageData,
			 this->GradientImageIsGradient,
                         this->OpacityTF, this->ColourTF,
                         this->OmegaL, this->OmegaH, min_gradmag,
                         &vectorDx, this->Px,
                         &vectorDy, this->Py,
                         &vectorDz, this->Pz);
        }
        break;
        default:
            vtkErrorMacro(<< "Unknown type of scalar data!");
        }

        // now we should copy vectorD over to D
        if (this->Dx)
            delete this->Dx;

        this->Dx = new ShellVoxel[vectorDx.size()];
        for (unsigned i = 0; i < vectorDx.size(); i++)
            memcpy(this->Dx + i, &(vectorDx[i]), sizeof(ShellVoxel));

        // now we should copy vectorD over to D
        if (this->Dy)
            delete this->Dy;

        this->Dy = new ShellVoxel[vectorDy.size()];
        for (unsigned i = 0; i < vectorDy.size(); i++)
            memcpy(this->Dy + i, &(vectorDy[i]), sizeof(ShellVoxel));

        // now we should copy vectorD over to D
        if (this->Dz)
            delete this->Dz;

        this->Dz = new ShellVoxel[vectorDz.size()];
        for (unsigned i = 0; i < vectorDz.size(); i++)
            memcpy(this->Dz + i, &(vectorDz[i]), sizeof(ShellVoxel));
        
        cout << vectorDx.size() << " shell voxels found." << endl;

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

