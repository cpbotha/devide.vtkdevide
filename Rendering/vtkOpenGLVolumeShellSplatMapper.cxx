// vtkOpenGLVolumeShellSplatMapper copyright (c) 2003 
// by Charl P. Botha cpbotha@ieee.org 
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id: vtkOpenGLVolumeShellSplatMapper.cxx,v 1.20 2004/01/08 15:00:52 cpbotha Exp $
// vtk class for volume rendering by shell splatting

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

#if defined(WIN32)
// it seems Visual C++ math.h doesn't declare this (schmuckheads)
   #ifndef M_PI
      #define M_PI 3.14159265358979323846
   #endif
// and due to Visual C++ 6.0 "language extensions" (actually M$ raping the
// standard: http://support.microsoft.com/default.aspx?scid=KB;en-us;q167748)
// we need this workaround:
   #define for if(0);else for
#endif

#include <time.h>
#include "vtkCamera.h"
#include "vtkColorTransferFunction.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPiecewiseFunction.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkOpenGLVolumeShellSplatMapper.h"
#include "vtkRenderer.h"
#include "vtkShellExtractor.h"
#include "vtkTransform.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#ifdef __APPLE__
   #include <OpenGL/gl.h>
#else
   #include <GL/gl.h>
#endif


// discrete Gaussian will be N x N
#define OGLVSM_RF_N 64
#define SSM_VERBOSE_OUTPUT

// stupid global variable we can use
long voxels_drawn;

void (*DrawVoxel)(ShellVoxel*,
                  const unsigned char&, const int&, const int&,
                  GLfloat*,
                  const GLfloat&, const GLfloat&,
                  GLfloat*, GLfloat*);

void DrawVoxelSplat(ShellVoxel* Dptr,
                    const unsigned char& octantIdx, const int& y, const int& z,
                    GLfloat* prev_colour,
                    const GLfloat& ambient, const GLfloat& diffuse,
                    GLfloat* u, GLfloat* v)
{
   static GLfloat temp_mat[4];
//   static GLfloat vt_array[12];
   static GLfloat ambiento;
   static GLfloat diffuseo;
   if (vtkShellExtractor::shell_noc_visibility_lut[octantIdx][Dptr->nbrOCode])
       //if (1)
   {
      if (Dptr->Red != prev_colour[0] || Dptr->Green != prev_colour[1] ||
          Dptr->Blue != prev_colour[2] || Dptr->Opacity != prev_colour[3])
      //if (0)
      {
         prev_colour[0] = Dptr->Red;
         prev_colour[1] = Dptr->Green;
         prev_colour[2] = Dptr->Blue;
         prev_colour[3] = Dptr->Opacity;
         // we premultiply colour with opacity, because that's what Porter and Duff say
         ambiento = ambient * Dptr->Opacity;
         temp_mat[0] = prev_colour[0] * ambiento;
         temp_mat[1] = prev_colour[1] * ambiento;
         temp_mat[2] = prev_colour[2] * ambiento;
         temp_mat[3] = prev_colour[3];
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, temp_mat);
         diffuseo = diffuse * Dptr->Opacity;
         temp_mat[0] = prev_colour[0] * diffuseo;
         temp_mat[1] = prev_colour[1] * diffuseo;
         temp_mat[2] = prev_colour[2] * diffuseo;
         // we don't need to set temp_mat[3] again, it's already done
         // the diffuse opacity is the resultant lit opacity
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, temp_mat);
      }

      glNormal3fv(Dptr->Normal);


      // change to:
      // glInterleavedArrays(GL_T2F_V3F, stride=0, pointer), i.e.
      // 2 float texture coords followed by 3 float vertices and
      // glDrawArrays(GL_QUADS, 0, 4)

      /*
       vt_array[0] = Dptr->x - u[0] - v[0];
       vt_array[1] = y - u[1] - v[1];
       vt_array[2] = z - u[2] - v[2];
       * 
       vt_array[3] = Dptr->x + u[0] - v[0]; 
       vt_array[4] = y + u[1] - v[1]; 
       vt_array[5] = z + u[2] - v[2];
       * 
       vt_array[6] = Dptr->x + u[0] + v[0]; 
       vt_array[7] = y + u[1] + v[1]; 
       vt_array[8] = z + u[2] + v[2];
       * 
       vt_array[9] = Dptr->x - u[0] + v[0]; 
       vt_array[10] = y - u[1] + v[1]; 
       vt_array[11] = z - u[2] + v[2];
       * 
       glInterleavedArrays(GL_V3F, 0, vt_array);
       glDrawArrays(GL_QUADS, 0, 4);
       */

      glTexCoord2f(0.0, 0.0);
      glVertex3f(Dptr->x - u[0] - v[0], y - u[1] - v[1], z - u[2] - v[2]);

      glTexCoord2f(1.0, 0.0);
      glVertex3f(Dptr->x + u[0] - v[0], y + u[1] - v[1], z + u[2] - v[2]);

      glTexCoord2f(1.0, 1.0);
      glVertex3f(Dptr->x + u[0] + v[0], y + u[1] + v[1], z + u[2] + v[2]);

      glTexCoord2f(0.0, 1.0);
      glVertex3f(Dptr->x - u[0] + v[0], y - u[1] + v[1], z - u[2] + v[2]);

      voxels_drawn++;
   }
}


void DrawVoxelPoint(ShellVoxel* Dptr,
                    const unsigned char& octantIdx, const int& y, const int& z,
                    GLfloat* prev_colour,
                    const GLfloat& ambient, const GLfloat& diffuse,
                    GLfloat* u, GLfloat* v)
{
   static GLfloat temp_mat[4];
   if (vtkShellExtractor::shell_noc_visibility_lut[octantIdx][Dptr->nbrOCode])
   {
      if (Dptr->Red != prev_colour[0] || Dptr->Green != prev_colour[1] ||
          Dptr->Blue != prev_colour[2] || Dptr->Opacity != prev_colour[3])
      {
         prev_colour[0] = Dptr->Red;
         prev_colour[1] = Dptr->Green;
         prev_colour[2] = Dptr->Blue;
         prev_colour[3] = Dptr->Opacity;
         temp_mat[0] = prev_colour[0] * ambient;
         temp_mat[1] = prev_colour[1] * ambient;
         temp_mat[2] = prev_colour[2] * ambient;
         temp_mat[3] = prev_colour[3];
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, temp_mat);
         temp_mat[0] = prev_colour[0] * diffuse;
         temp_mat[1] = prev_colour[1] * diffuse;
         temp_mat[2] = prev_colour[2] * diffuse;
         // we don't need to set temp_mat[3] again, it's already done
         // the diffuse opacity is the resultant lit opacity
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, temp_mat);
      }

      glNormal3fv(Dptr->Normal);
      glVertex3f(Dptr->x, y, z);

      voxels_drawn++;
   }
}

// sigma == standard deviation
void IntegratedGaussian(double* output, int len, float cutoff, float sigma)
{
   // calculate sample step
   float inc = cutoff * 2 / (float)len;
   // pre-calculate 1 / (sqrt(2 * M_PI) * sigma^3)
   float mfactor = 1.0  / (2.0 * M_PI * sigma * sigma);
   float x, y;
   float xsqr, ysqr;
   float cutoffsqr = cutoff * cutoff;
   int idx = 0;


   for (int iy = 0; iy < len; iy++)
   {
      y = -cutoff + (float)iy * inc;
      ysqr = y * y;
      for (int ix = 0; ix < len; ix++)
      {
         x = -cutoff + (float)ix * inc;
         xsqr = x*x;

         if (xsqr + ysqr < cutoffsqr)
            //if (1)
            output[idx] = mfactor * exp(- (xsqr + ysqr) /
                                        (2.0 * sigma * sigma));
         else
            output[idx] = 0.0;

         // KLUDGE ALERT, FIXME -- whoop, no kludge anymore (I think)
         //output[idx] /= mfactor;

         idx++;
      }
   }

}

vtkStandardNewMacro(vtkOpenGLVolumeShellSplatMapper);

vtkOpenGLVolumeShellSplatMapper::vtkOpenGLVolumeShellSplatMapper()
{
   //this->Input = NULL;
   this->ShellExtractor = vtkShellExtractor::New();

   this->volMatrix = vtkMatrix4x4::New();
   this->VoxelsToWorldMatrix = vtkMatrix4x4::New();
   this->WorldToVoxelsMatrix = vtkMatrix4x4::New();
   this->ViewMatrix = vtkMatrix4x4::New();
   this->VoxelsToViewMatrix = vtkMatrix4x4::New();
   this->ViewToVoxelsMatrix = vtkMatrix4x4::New();
   this->PerspectiveMatrix = vtkMatrix4x4::New();

   this->integrated_rfunc = new double[OGLVSM_RF_N * OGLVSM_RF_N];
   // for each and every view, integrated_rfunc will be normalised into this - and we will add
   // alpha values to make this a luminance-alpha buffer
   this->normalised_integrated_rfunc = new float[OGLVSM_RF_N * OGLVSM_RF_N * 2];

   this->GradientImageData = NULL;

   // this will set the Gaussian parameters and initialise the textures themselves
   this->SetRenderMode(0);

   // this will setup the initial rendering mode (to normal PBTF)
   this->SetPerspectiveOrderingMode(0);
}

vtkOpenGLVolumeShellSplatMapper::~vtkOpenGLVolumeShellSplatMapper()
{
   delete this->integrated_rfunc;
   delete this->normalised_integrated_rfunc;

   // kill a bunch of matrices and transforms
   volMatrix->Delete();
   VoxelsToWorldMatrix->Delete();
   WorldToVoxelsMatrix->Delete();
   ViewMatrix->Delete();
   VoxelsToViewMatrix->Delete();
   ViewToVoxelsMatrix->Delete();
   PerspectiveMatrix->Delete();

   this->SetGradientImageData(NULL);
   this->SetInput((vtkImageData*)NULL);
   if (ShellExtractor)
      ShellExtractor->Delete();
}

void vtkOpenGLVolumeShellSplatMapper::CalculatePerViewMatrices(vtkRenderer* ren, vtkVolume* vol)
{
   // the vtkVolume matrix is a volume to world matrix, get it
   this->volMatrix->DeepCopy(vol->GetMatrix());

   // transform to account for origin and spacing of input volume
   vtkTransform *voxelsTransform = vtkTransform::New();
   voxelsTransform->Identity();
   voxelsTransform->Translate(this->GetInput()->GetOrigin());
   voxelsTransform->Scale(this->GetInput()->GetSpacing());
   // default is premultiply, so now we have:
   // translate * scale * [point to transform from vox to volume]

   // volumetoworldmatrix * voxeltovolumematrix == voxeltoworldmatrix
   vtkTransform *voxelsToViewTransform = vtkTransform::New();
   voxelsToViewTransform->SetMatrix(this->volMatrix);
   // NOTE: premultiply() means that what we we add from now has to be
   // applied to the transformee (not available yet) FIRST, i.e. 
   // new_transformation_matrix = previous_transformation_matrix * added_operation
   voxelsToViewTransform->PreMultiply();
   voxelsToViewTransform->Concatenate(voxelsTransform->GetMatrix());
   // so, now we have: volumetoworld * translate * scale [* some_point]
   // i.e. at this point it's more a voxels to world transform

   // then we need a voxelstoworldmatrix
   // volumetoworldmatrix * translate * scale == voxeltoworldmatrix
   this->VoxelsToWorldMatrix->DeepCopy(voxelsToViewTransform->GetMatrix());

   // copy the voxeltoworldmatrix
   this->WorldToVoxelsMatrix->DeepCopy(voxelsToViewTransform->GetMatrix());
   // invert to get worldtovoxelmatrix
   this->WorldToVoxelsMatrix->Invert();

   // now continue with the voxelstoviewtransform
   this->ViewMatrix->DeepCopy(ren->GetActiveCamera()->GetViewTransformMatrix());

   // this was volumetoworld * translate * scale [* some_voxel_point]
   voxelsToViewTransform->PostMultiply();
   // now it's:
   // worldtoview * volumetoworld * translate * scale [* some_voxel_point]
   voxelsToViewTransform->Concatenate(this->ViewMatrix);
   this->VoxelsToViewMatrix->DeepCopy(voxelsToViewTransform->GetMatrix());

   this->ViewToVoxelsMatrix->DeepCopy(this->VoxelsToViewMatrix);
   this->ViewToVoxelsMatrix->Invert();

   // get perspectivematrix from current camera (also using aspect from renderer)
   ren->ComputeAspect();
   //this->PerspectiveMatrix->DeepCopy(ren->GetActiveCamera()->GetPerspectiveTransformMatrix(ren->GetAspect()[0] / ren->GetAspect()[1], -0.1, 0.1));
   // we don't want to bring the aspect into this, else our quad construction goes tits-up (I'm not entirely sure why)
   this->PerspectiveMatrix->DeepCopy(ren->GetActiveCamera()->GetPerspectiveTransformMatrix(1.0, -0.1, 0.1));

   voxelsTransform->Delete();
   voxelsToViewTransform->Delete();

}

void vtkOpenGLVolumeShellSplatMapper::SetGaussianRadialExtent(double gre)
{
   this->gaussian_radial_extent = gre;
   IntegratedGaussian(this->integrated_rfunc, OGLVSM_RF_N, 
                      this->gaussian_radial_extent, this->gaussian_sigma); // sigma was 0.7
}

void vtkOpenGLVolumeShellSplatMapper::SetGaussianSigma(double gs)
{
   this->gaussian_sigma = gs;
   IntegratedGaussian(this->integrated_rfunc, OGLVSM_RF_N, 
                      this->gaussian_radial_extent, this->gaussian_sigma); // sigma was 0.7
}

void vtkOpenGLVolumeShellSplatMapper::SetRenderMode(int newRenderMode)
{
   this->RenderMode = newRenderMode;

   if (RenderMode == 0)
   {
      this->EllipsoidDiameter = 4.0;
      this->gaussian_radial_extent = 2.0;
      this->gaussian_sigma = 0.7;
   }
   else
   {
      // in modes 1 and 2 we're not texturing, so Gaussian params aren't important
      this->EllipsoidDiameter = 1.6;
   }

   IntegratedGaussian(this->integrated_rfunc, OGLVSM_RF_N, 
                      this->gaussian_radial_extent, this->gaussian_sigma); // sigma was 0.7

}

void vtkOpenGLVolumeShellSplatMapper::Render(vtkRenderer* ren, vtkVolume* vol)
{
   if (this->GetInput() == NULL)
   {
      vtkErrorMacro(<< "No input set!");
      return;
   }
   else
   {
      this->GetInput()->UpdateInformation();
      this->GetInput()->SetUpdateExtentToWholeExtent();
      this->GetInput()->Update();
   }

   //   if (!ren->GetActiveCamera()->GetParallelProjection())
   //   {
   //      vtkWarningMacro(<< "ShellSplatter does not support perspective rendering yet.  Please do e.g.: ren->GetActiveCamera()->ParallelProjectionOn().");
   //   }

#ifdef SSM_VERBOSE_OUTPUT   
   clock_t start_clock = clock();
#endif

   // these can't be NULL, since vtkVolumeProperty will MAKE them if
   // they don't exist
   vtkColorTransferFunction* ColourTF = 
   vol->GetProperty()->GetRGBTransferFunction(0);
   vtkPiecewiseFunction* OpacityTF =
   vol->GetProperty()->GetScalarOpacity(0);

   // all the Set[Object]Macros first check whether the passed
   // parameter differs from the variable, so we don't have to
   // worry that the ShellExtractor is going to run unnecessarily
   ShellExtractor->SetInput(this->GetInput());
   ShellExtractor->SetGradientImageData(this->GradientImageData);
   ShellExtractor->SetOpacityTF(OpacityTF);
   ShellExtractor->SetColourTF(ColourTF);
   ShellExtractor->SetOmegaL(this->OmegaL);
   ShellExtractor->SetOmegaH(this->OmegaH);
   // this should only happen when necessary
   ShellExtractor->Update();

   CalculatePerViewMatrices(ren, vol);


   // get D and P in temp pointers
   ShellVoxel* D = (ShellVoxel*)(this->ShellExtractor->GetD());
   int* P = this->ShellExtractor->GetP();


   // woef, let's render
   // first make sure we're working in the right context
   ren->GetRenderWindow()->MakeCurrent();

   // store all gl attributes for this context
   glPushAttrib(GL_ALL_ATTRIB_BITS);

   // first activate blending
   glEnable(GL_BLEND);

   if (this->RenderMode == 2)
   {
      // point rendering does not have premultiplied alpha
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }
   else
   {
      // mapped and non-mapped quad
      glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
   }

   // so we can blend with geometry - VTK does opaque -> translucent -> overlay!
   glEnable(GL_DEPTH_TEST); // if there's an opaque object closer, we can not paint there
   glDepthMask(GL_FALSE); // but our thingies do NOT affect the z-buffer
   // make sure lighting is on
   glEnable(GL_LIGHTING);
   // we want a single shaded colour per polygon
   glShadeModel(GL_FLAT);
   // we're not interested in two-sided lighting, and I hope this makes it 
   // a bit faster) - it also solves complications with our quad not always
   // having the same vertex order
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

   GLfloat ambient = vol->GetProperty()->GetAmbient();
   GLfloat diffuse = vol->GetProperty()->GetDiffuse();

   GLfloat temp_mat[4];
   temp_mat[3] = 1.0;

   GLfloat specular = vol->GetProperty()->GetSpecular();
   for (int i = 0; i < 3; i++)
   {
      temp_mat[i] = specular;   
   }
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, temp_mat);

   GLfloat shininess = vol->GetProperty()->GetSpecularPower();
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

   vtkMatrix4x4 *Q = vtkMatrix4x4::New();
   Q->Zero();
   // to get quadric matrix diagonals: 1/((desired_diameter/2)^2)
   // we want to create a quadric that is centered at the origin, and has
   // a diameter (e.g. 2a) equal to OGLVSM_VOXEL_DIAMETER
   double half_diameter = this->EllipsoidDiameter / 2.0; // calc a
   double temp_elem = 1.0 / (half_diameter * half_diameter); // then 1 / a^2
   Q->SetElement(0,0,temp_elem);
   Q->SetElement(1,1,temp_elem);
   Q->SetElement(2,2,temp_elem);
   Q->SetElement(3,3,-1);

   // Qa = inv(M)' Q inv(M)
   // i.e. to transform Q with M, we multiply inverse of M transposed with Q with inverse of M
   vtkMatrix4x4 *vtvmt = vtkMatrix4x4::New();
   vtkMatrix4x4 *tempm = vtkMatrix4x4::New();
   vtkMatrix4x4 *Qa = vtkMatrix4x4::New();
   vtkMatrix4x4 *M = vtkMatrix4x4::New();
   vtkMatrix4x4 *Mi = vtkMatrix4x4::New();

   //M->DeepCopy(this->VoxelsToViewMatrix);
   vtkMatrix4x4::Multiply4x4(this->PerspectiveMatrix, this->VoxelsToViewMatrix,  M);
   // no translation thank you - yes we can remove this as we put it there
   // in the first place
   M->Element[0][3] = 0.0;
   M->Element[1][3] = 0.0;
   M->Element[2][3] = 0.0;
   // FIXME: we want a normalised camera, and this might be stopping us
   // we could try the full: cA - bb' instead of just killing translation

   vtkDebugMacro( << *M << endl );
   vtkMatrix4x4::Invert(M, Mi);

   vtkMatrix4x4::Multiply4x4(Q, Mi, tempm);
   vtkMatrix4x4::Transpose(Mi, vtvmt);
   vtkMatrix4x4::Multiply4x4(vtvmt, tempm, Qa);

   vtkDebugMacro( << *Qa << endl );

   // calculate b * b' yielding the 3x3 bb'
   /*
      cout << "bb == " << endl;
      double bb[3][3];
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j <3; j++)
         {
            bb[i][j] = Qa->Element[i][3] * Qa->Element[3][j];
            Qa->Element[i][j] = Qa->Element[3][3] * Qa->Element[i][j] - bb[i][j];
            cout << bb[i][j] << " ";
         }
         cout << endl;
      }
   
      cout << *Qa << endl;
   */
   double *pQa[3];
   double Qa_eigvec[3][3];
   double *pQa_eigvec[3];
   double Qa_eigval[3];

   for (int i = 0; i < 3; i++)
   {
      pQa[i] = Qa->Element[i];
      pQa_eigvec[i] = Qa_eigvec[i];
   }

   vtkMath::JacobiN(pQa, 3, Qa_eigval, pQa_eigvec);

   //   cout << "Qa_eigval == ";
   //   for (int i = 0; i < 3; i++)
   //   {
   //      cout << Qa_eigval[i] << " (";
   //      for (int j = 0; j < 3; j++)
   //      {
   //         cout << Qa_eigvec[j][i] << " ";
   //      }
   //      cout << ") ";
   //   }
   //   cout << endl << endl;

   // VERY IMPORTANT:
   // we only use the FIRST THREE eigvecs, BUT:
   // now we have to arrange the matrices so that the first two eigvecs 
   // have 0 as their third element and the last has 1 as its last element
   // see: "Model-based 3d tracking of an articulated hand" by 
   // Stenger, Mendonca, Cipolla
   // to make a robust implementation, we check for the eigvec with the
   // largest third element (i.e. approximately 1.0)

   GLfloat u[4], v[4];

   int largest_idx = 0;
   double largest_third = 0.0;

   for (int i = 0; i < 3; i++)
   {
      if (fabs(Qa_eigvec[2][i]) > largest_third)
      {
         largest_third = fabs(Qa_eigvec[2][i]);
         largest_idx = i;
      }
      u[i] = v[i] = 0.0;
   }

   int used_idxs[2];
   for (int i = 0, j = 0; i < 3; i++)
   {
      if (i != largest_idx)
      {
         used_idxs[j] = i;
         j++;
      }
   }

   // the 0.5 is here for clarity: we've transformed a quadric of the
   // same dimensions as a voxel, but now we want a u and a v that
   // are going to be added twice (e.g. p - u -v to p + u + v) in the
   // quad rendering procedure
   double mag0 = 0.5 * 2.0 * sqrt(1 / fabs(Qa_eigval[used_idxs[0]]));
   double mag1 = 0.5 * 2.0 * sqrt(1 / fabs(Qa_eigval[used_idxs[1]]));
   for (int i = 0; i < 2; i++)
   {
      u[i] =  mag0 * Qa_eigvec[i][used_idxs[0]];
      v[i] =  mag1 * Qa_eigvec[i][used_idxs[1]];
   }


   vtkDebugMacro( << "u = " << u[0] << " " << u[1] << " " << u[2] << " " << u[3]);
   vtkDebugMacro( << "v = " << v[0] << " " << v[1] << " " << v[2] << " " << v[3]);

   Mi->MultiplyPoint(u,u);
   Mi->MultiplyPoint(v,v);


   // let's normalise the Gaussian
   double gsum = 0.0;
   for (int i = 0; i < OGLVSM_RF_N * OGLVSM_RF_N; i++)
   {
      gsum += this->integrated_rfunc[i];
   }

   // now transform u and v to world coordinates from voxel coords
   vtkMatrix4x4 *v2wnt = vtkMatrix4x4::New();
   v2wnt->DeepCopy(this->VoxelsToWorldMatrix);
   v2wnt->Element[0][3] = v2wnt->Element[1][3] = v2wnt->Element[2][3] = 0.0;
   float uw[4], uwmag, vw[4], vwmag;
   v2wnt->MultiplyPoint(u, uw);
   uwmag = vtkMath::Norm(uw);
   v2wnt->MultiplyPoint(v, vw);
   vwmag = vtkMath::Norm(vw);
   v2wnt->Delete();

   // this one should work! this is if the Gaussian is created with mfactor (i.e. NO kludge) - we scale our values by how 
   // much our circular area (pi * r^2) differs from the "total" area (i.e. area of circle with radius == 4 * sigma, where 
   // the Gaussian is normally 0.0)
   double nfactor = (this->gaussian_radial_extent * this->gaussian_radial_extent) / (4.0 * this->gaussian_sigma * 4.0 * this->gaussian_sigma);
#ifdef SSM_VERBOSE_OUTPUT
   cout << "nfactor == " << nfactor << endl;
#endif
   for (int i = 0; i < OGLVSM_RF_N * OGLVSM_RF_N; i++)
   {
      // uncomment below to have both opacity and  colour modulated by the Gaussian
      this->normalised_integrated_rfunc[i * 2] = this->normalised_integrated_rfunc[i * 2 + 1] = (float)(this->integrated_rfunc[i] / nfactor);
      // uncomment below to have only opacity modulated by the Gaussian - one
      // argue that the energy of a voxel is determined by its opacity
      //this->normalised_integrated_rfunc[i * 2] = 1.0;
      //this->normalised_integrated_rfunc[i * 2 + 1] = (float)(this->integrated_rfunc[i] / nfactor);

   }


   float u2mag = 0.0;

   for (int i = 0; i < 4; i++)
   {
      u2mag += u[i] * u[i];
   }

   vtkDebugMacro( << "u2 = " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << " mag == " << sqrt(u2mag));

   float v2mag = 0.0;
   for (int i = 0; i < 4; i++)
   {
      v2mag += v[i] * v[i];
   }

   vtkDebugMacro( << "v2 = " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " mag == " << sqrt(v2mag));

   Q->Delete();
   vtvmt->Delete();
   tempm->Delete();
   Qa->Delete();
   M->Delete();
   Mi->Delete();

   // let's try rendering in voxel world
   double *mat = VoxelsToViewMatrix->Element[0];

   double mat2[16];
   mat2[0] = mat[0];
   mat2[1] = mat[4];
   mat2[2] = mat[8];
   mat2[3] = mat[12];
   mat2[4] = mat[1];
   mat2[5] = mat[5];
   mat2[6] = mat[9];
   mat2[7] = mat[13];
   mat2[8] = mat[2];
   mat2[9] = mat[6];
   mat2[10] = mat[10];
   mat2[11] = mat[14];
   mat2[12] = mat[3];
   mat2[13] = mat[7];
   mat2[14] = mat[11];
   mat2[15] = mat[15];
   // insert model transformation, now we can draw in voxel coords
   // we've taken over the complete model view matrix... the projection
   // matrix we're leaving for now (to the camera); usually the camera does
   // world -> view as well
   glMatrixMode( GL_MODELVIEW );
   glPushMatrix();
   glLoadIdentity();
   glMultMatrixd(mat2);

   long tex_idx = -1;

   if (RenderMode == 0)
   {

#ifdef GL_VERSION_1_1
     // we have to do this to unbind whatever a vtkTexture has bound to
     // this target (GL_TEXTURE_2D), although we're never going to rebind it;
     // we delete it at the end of this method
     // if we didn't bind a new texture to this target, our uploads would
     // overwrite whichever texture was bound to that target when we started
     // mucking around... remember that these newfangled texture objects are
     // mutable
     GLuint tempIndex = 0;
     glGenTextures(1, &tempIndex);
     tex_idx = (long) tempIndex;
     glBindTexture(GL_TEXTURE_2D, tempIndex);
#else
     tex_idx = glGenLists(1);
     glDeleteLists ((GLuint) tex_idx, (GLsizei) 0);
     glNewList ((GLuint) tex_idx, GL_COMPILE);
#endif

     glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                       GL_NEAREST);
     glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                       GL_NEAREST );

      // if internal format is  GL_INTENSITY: C = Cf * I, A = Af * I
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


      // modulate colour AND opacity with the gaussian
      glTexImage2D( GL_TEXTURE_2D, 0 , 4, // 4 components in texture
                    OGLVSM_RF_N, OGLVSM_RF_N, 0, GL_LUMINANCE_ALPHA, // format val0 -> r,g,b and val1 -> opacity
                    GL_FLOAT, (const GLvoid *) this->normalised_integrated_rfunc );

#ifdef GL_VERSION_1_1
      // do we really need this call, i.e. is the texture not already
      // bound?
      glBindTexture(GL_TEXTURE_2D, tex_idx);
#else
      glEndList();
      glCallList ((GLuint) tex_idx);
#endif


      glEnable(GL_TEXTURE_2D);

      glBegin(GL_QUADS);
      DrawVoxel = DrawVoxelSplat;
   }
   else if (this->RenderMode == 1)
   {
      // in this rendermode, we don't do any texture mapping, only quads
      glBegin(GL_QUADS);
      DrawVoxel = DrawVoxelSplat;
   }
   else
   {
      float largest_mag = mag0 > mag1 ? mag0 : mag1;


      float max_siz = largest_mag * ren->GetSize()[0] / 2.0;
      if (max_siz < 1.0)
         max_siz = 1.0;

      //glPointSize(2.0 * max_siz);
      //glEnable(GL_POINT_SMOOTH);
      glBegin(GL_POINTS);
      DrawVoxel = DrawVoxelPoint;
   }


   // we need this for the iteration!
   int inputDims[3];
   int  xdim, ydim, zdim;
   this->GetInput()->GetDimensions(inputDims);
   xdim = inputDims[0];
   ydim = inputDims[1];
   zdim = inputDims[2];

   // we have to find out WHERE the view point is (once)
   double camWorldPos[4];
   ren->GetActiveCamera()->GetPosition(camWorldPos);
#ifdef SSM_VERBOSE_OUTPUT
   cout << "camWorldPos == " << camWorldPos[0] << "," << camWorldPos[1] << "," << camWorldPos[2] << endl;
#endif
   camWorldPos[3] = 1.0;
   // and convert it to voxel space
   double camVoxelPos[4];
   this->WorldToVoxelsMatrix->MultiplyPoint(camWorldPos, camVoxelPos);

   int octantIdx;
   voxels_drawn = 0;

   if (ren->GetActiveCamera()->GetParallelProjection())
   {
      // special case if this is an orthogonal projection:
      // we need to project the centre of the volume onto the view plane and
      // then check in which volume octant that point is

      vtkMatrix4x4 *octantM = vtkMatrix4x4::New();
      vtkMatrix4x4::Multiply4x4(this->PerspectiveMatrix, this->VoxelsToViewMatrix,  octantM);
      float voxel_volume_centre[4], projected_volume_centre[4], projected_voxel_volume_centre[4];
      voxel_volume_centre[0] = (float)xdim / 2.0;
      voxel_volume_centre[1] = (float)ydim / 2.0;
      voxel_volume_centre[2] = (float)zdim / 2.0;
      voxel_volume_centre[3] = 1.0;

      // project the volume centre onto the view plane
      octantM->MultiplyPoint(voxel_volume_centre, projected_volume_centre);

      // set its third dimension to zero (that should be IN the plane)
      projected_volume_centre[2] = 0.0;

      // we're going back to voxel space
      octantM->Invert();
      octantM->MultiplyPoint(projected_volume_centre, projected_voxel_volume_centre);

      // now determine the octant
      // and yes, it does make sense that we're checking with <: the projected
      // point is not stuck on the view-plane, it's slipping and sliding as
      // the view plane is moving.  Think on that and take it from there. :)
      octantIdx = 0;
      if (projected_voxel_volume_centre[0] < voxel_volume_centre[0])
      {
         octantIdx |= 0x1;
      }
      if (projected_voxel_volume_centre[1] < voxel_volume_centre[1])
      {
         octantIdx |= 0x2;
      }
      if (projected_voxel_volume_centre[2] < voxel_volume_centre[2])
      {
         octantIdx |= 0x4;
      }
      
#ifdef SSM_VERBOSE_OUTPUT
      cout << "PVC " << projected_voxel_volume_centre[0] << " " << projected_voxel_volume_centre[1] << " " << projected_voxel_volume_centre[2] << " " << endl;
#endif

      // take care of the matrix we made
      octantM->Delete();

      this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2], octantIdx, D, P, inputDims[1], ambient, diffuse, u, v);
   }
   else if (this->PerspectiveOrderingMode == 0) // PERSPECTIVE rendering, PBTF
   {
      // check whether it's volume-on, face-on, edge-on or corner-on
      int xin = 0, yin = 0, zin = 0, xyztot = 0;
      if (camVoxelPos[0] >= 0 && camVoxelPos[0] < inputDims[0])
      {
         xin = 1;
         xyztot++;
      }
      if (camVoxelPos[1] >= 0 && camVoxelPos[1] < inputDims[1])
      {
         yin = 1;
         xyztot++;
      }
      if (camVoxelPos[2] >= 0 && camVoxelPos[2] < inputDims[2])
      {
         zin = 1;
         xyztot++;
      }

      int pq, pr, ps;

      // **************************************************************************
      // NB: convention: LDB is our origin, NOT LDF as in Swan's thesis!
      // **************************************************************************

      // volume-on: break up volume into 8 octants
      if (xyztot == 3)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "VOLUME-ON" << endl;
#endif         
         pq = vtkMath::Round(camVoxelPos[0]);
         pq = (pq >= xdim) ? xdim - 1 : pq;
         pr = vtkMath::Round(camVoxelPos[1]);
         pr = (pr >= ydim) ? ydim - 1 : pr;
         ps = vtkMath::Round(camVoxelPos[2]);
         ps = (ps >= zdim) ? zdim - 1 : ps;
         // rendering octant LDF (code: 4)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      ps + 1, zdim,   3, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDF (code: 5)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      ps + 1, zdim,   2, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUF (code: 6)      
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   ps + 1, zdim,   1, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUF (code: 7)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   ps + 1, zdim,   0, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LDB (code: 0)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      0, ps + 1,      7, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDB (code: 1)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      0, ps + 1,      6, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUB (code: 2)
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   0, ps + 1,      5, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUB (code: 3)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   0, ps + 1,      4, D, P, ydim, ambient, diffuse, u, v);
      } // end VOLUME-ON

      // face-on: break volume into 4 quadrants
      else if (xyztot == 2)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "FACE-ON: ";
#endif         
         if (zin == 0)
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            if (camVoxelPos[2] >= zdim)
            {
#ifdef SSM_VERBOSE_OUTPUT                
               cout << "z >= zdim" <<endl;
#endif               
               // quadrant LU
               this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   0, zdim,   5, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, pq + 1,      0, pr + 1,      0, zdim,   7, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   0, zdim,   4, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      0, zdim,   6, D, P, inputDims[1], ambient, diffuse, u, v);
            }
            else
            {
#ifdef SSM_VERBOSE_OUTPUT                
               cout << "z < 0" <<endl;
#endif               
               // MIRROR of clause above, only shell rendering octantIdx changes
               // quadrant LU
               this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   0, zdim,   1, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, pq + 1,      0, pr + 1,      0, zdim,   3, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   0, zdim,   0, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      0, zdim,   2, D, P, inputDims[1], ambient, diffuse, u, v);
            }
         }
         else if (yin == 0)
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;
            if (camVoxelPos[1] >= ydim)
            {
               // We're viewing from "ABOVE"
               // quadrant LU
               this->DrawVoxels(0, pq + 1,        0, ydim,   0, ps + 1,        7, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, pq + 1,        0, ydim,   ps + 1, zdim,   3, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(pq + 1, xdim,   0, ydim,   0, ps + 1,        6, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(pq + 1, xdim,   0, ydim,   ps + 1, zdim,   2, D, P, inputDims[1], ambient, diffuse, u, v);
            }
            else
            {
               // MIRROR
               // quadrant LU
               this->DrawVoxels(0, pq + 1,        0, ydim,   0, ps + 1,        5, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, pq + 1,        0, ydim,   ps + 1, zdim,   1, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(pq + 1, xdim,   0, ydim,   0, ps + 1,        4, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(pq + 1, xdim,   0, ydim,   ps + 1, zdim,   0, D, P, inputDims[1], ambient, diffuse, u, v);
            }
         }
         else // xin == 0
         {
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;
            if (camVoxelPos[0] >= xdim)
            {
               // We're viewing from the "RIGHT"
               // quadrant LU
               this->DrawVoxels(0, xdim,   pr + 1, ydim,   ps + 1, zdim,   1, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, xdim,   0, pr + 1,      ps + 1, zdim,   3, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(0, xdim,   pr + 1, ydim,   0, ps + 1,      5, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(0, xdim,   0, pr + 1,      0, ps + 1,      7, D, P, inputDims[1], ambient, diffuse, u, v);
            }
            else
            {
               // MIRROR
               // quadrant LU
               this->DrawVoxels(0, xdim,   pr + 1, ydim,   ps + 1, zdim,   0, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant LD
               this->DrawVoxels(0, xdim,   0, pr + 1,      ps + 1, zdim,   2, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RU
               this->DrawVoxels(0, xdim,   pr + 1, ydim,   0, ps + 1,      4, D, P, inputDims[1], ambient, diffuse, u, v);
               // quadrant RD
               this->DrawVoxels(0, xdim,   0, pr + 1,      0, ps + 1,      6, D, P, inputDims[1], ambient, diffuse, u, v);
            }
         }
      } // FACE-ON

      // edge-on
      else if (xyztot == 1)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "EDGE-ON" << endl;
#endif         
         int edge_idx = 0;
         if (zin)
         {
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;
            if (camVoxelPos[0] >= xdim)
               edge_idx |= 0x1;
            if (camVoxelPos[1] >= ydim)
               edge_idx |= 0x2;
            // this is sorted, really
            this->DrawVoxels(0, xdim,   0, ydim,   0, ps + 1,      edge_idx + 4, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(0, xdim,   0, ydim,   ps + 1, zdim,   edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
         else if (yin)
         {
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            if (camVoxelPos[0] >= xdim)
               edge_idx |= 0x1;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, xdim, 0, pr + 1, 0, zdim, edge_idx + 2, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(0, xdim, pr + 1, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
         else // xin
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            if (camVoxelPos[1] >= ydim)
               edge_idx |= 0x2;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, pq + 1,      0, ydim, 0, zdim, edge_idx + 1, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(pq + 1, xdim,   0, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
      } // EDGE-ON

      // corner-on
      else
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "CORNER-ON" << endl;
#endif         
         // being here means Vp is by definition outside the volume:
         octantIdx = 0;
         if (camVoxelPos[0] >= inputDims[0])
         {
            octantIdx |= 0x1;
         }
         if (camVoxelPos[1] >= inputDims[1])
         {
            octantIdx |= 0x2;
         }
         if (camVoxelPos[2] >= inputDims[2])
         {
            octantIdx |= 0x4;
         }
         this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2], octantIdx, D, P, inputDims[1], ambient, diffuse, u, v);
      } // CORNER-ON

   } // else PERSPECTIVE rendering, PBTF


   // ------------------------------------------------------------------------
   else if (this->PerspectiveOrderingMode == 1) // PERSPECTIVE rendering, iPBTF
   {
       // variables used for traversal of P and D
       int Pidx, Pidx0, Pidx1; // octant-numbering scheme
       ShellVoxel *Dptr, *Dptr0, *Dptr1;
       Dptr = Dptr0 = Dptr1 = (ShellVoxel*)NULL;

       // materials caching
       GLfloat prev_colour[4];
       prev_colour[0] = prev_colour[1] = prev_colour[2] = prev_colour[3] = -1;
       
      // check whether it's volume-on, face-on, edge-on or corner-on
      int xin = 0, yin = 0, zin = 0, xyztot = 0;
      if (camVoxelPos[0] >= 0 && camVoxelPos[0] < inputDims[0])
      {
         xin = 1;
         xyztot++;
      }
      if (camVoxelPos[1] >= 0 && camVoxelPos[1] < inputDims[1])
      {
         yin = 1;
         xyztot++;
      }
      if (camVoxelPos[2] >= 0 && camVoxelPos[2] < inputDims[2])
      {
         zin = 1;
         xyztot++;
      }

      int pq, pr, ps;

      // **************************************************************************
      // NB: convention: LDB is our origin, NOT LDF as in Swan's thesis!
      // **************************************************************************

      // volume-on: break up volume into 8 octants
      if (xyztot == 3)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "VOLUME-ON" << endl;
#endif         
         pq = vtkMath::Round(camVoxelPos[0]);
         pq = (pq >= xdim) ? xdim - 1 : pq;
         pr = vtkMath::Round(camVoxelPos[1]);
         pr = (pr >= ydim) ? ydim - 1 : pr;
         ps = vtkMath::Round(camVoxelPos[2]);
         ps = (ps >= zdim) ? zdim - 1 : ps;
         // rendering octant LDF (code: 4)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      ps + 1, zdim,   3, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDF (code: 5)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      ps + 1, zdim,   2, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUF (code: 6)      
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   ps + 1, zdim,   1, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUF (code: 7)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   ps + 1, zdim,   0, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LDB (code: 0)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      0, ps + 1,      7, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDB (code: 1)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      0, ps + 1,      6, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUB (code: 2)
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   0, ps + 1,      5, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUB (code: 3)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   0, ps + 1,      4, D, P, ydim, ambient, diffuse, u, v);
      } // end VOLUME-ON

      // face-on: break volume into 4 quadrants
      else if (xyztot == 2)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "FACE-ON: ";
#endif         
         if (zin == 0)
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;

            cout << "HINK: face-on, zin == 0, pq == "
                 << pq << " pr == " << pr << endl;

            
            // z is the slowest changing dimension
            int z0, z1, zinc;
            unsigned char octantIdxBYBX, octantIdxBYSX;
            unsigned char octantIdxSYBX, octantIdxSYSX;
            
            if (camVoxelPos[2] >= zdim)
            {
                z0 = 0;
                z1 = zdim;
                zinc = 1;
                // z is > zdim, so we can set this so long
                octantIdxBYBX =
                    octantIdxBYSX = octantIdxSYBX = octantIdxSYSX = 0x4;
            }
            else
            {
                z0 = zdim - 1;
                z1 = -1; // we're going to use != z1 as end cond
                zinc = -1;
                octantIdxBYBX =
                    octantIdxBYSX = octantIdxSYBX = octantIdxSYSX = 0;
            }

            // now partition y and x
            int bigyStart, bigyEnd, bigyInc, bigyThresh, smally, initsmally;

            if (pr >= (ydim - pr))
            {
                bigyStart = 0;
                bigyEnd = pr;
                bigyInc = 1;
                bigyThresh = pr - (ydim - pr); // biglen - small_len
                initsmally = ydim - 1;
                octantIdxBYBX |= 0x2; octantIdxBYSX = octantIdxBYBX;
            }
            else
            {
                bigyStart = ydim - 1;
                bigyEnd = pr - 1; // loop uses != bigyEnd
                bigyInc = -1;
                bigyThresh = pr * 2 - 1;
                initsmally = 0;
                octantIdxSYBX |= 0x2; octantIdxSYSX = octantIdxSYBX;
            }

            cout << bigyStart << " -> " << bigyEnd << " with "
                 << bigyInc << endl;
            cout << initsmally << " | " << bigyThresh << endl;

            
            int bigxStart, bigxEnd, bigxInc, bigxThresh, smallx, initsmallx;

            if (pq >= (xdim - pq))
            {
                bigxStart = 0;
                bigxEnd = pq;
                bigxInc = 1;
                bigxThresh = pq - (xdim - pq); // biglen - small_len
                initsmallx = xdim - 1;
                octantIdxBYBX |= 0x1; octantIdxSYBX = octantIdxBYBX;
            }
            else
            {
                bigxStart = xdim - 1;
                bigxEnd = pq - 1; // loop uses != bigxEnd
                bigxInc = -1;
                bigxThresh = pq * 2 - 1;
                initsmallx = 0;
                octantIdxBYSX |= 0x1; octantIdxSYSX = octantIdxBYSX;
            }

            
            bool yInterleaved;
            bool xInterleaved;

            ShellVoxel *DptrBY, *DptrBYBX, *DptrBYSX;
            DptrBY = DptrBYBX = DptrBYSX = (ShellVoxel*)NULL;
            
            ShellVoxel *DptrSY, *DptrSYBX, *DptrSYSX;
            DptrSY = DptrSYBX = DptrSYSX = (ShellVoxel*)NULL;
            
            int PidxBY, PidxSY;
            
            for (int z = z0; z != z1; z += zinc)
            {
                if (bigyStart == bigyThresh)
                    yInterleaved = true;
                else
                    yInterleaved = false;

                
                smally = initsmally;
                for (int bigy = bigyStart; bigy != bigyEnd; bigy += bigyInc)
                {
//                     cout << "y: " << bigy;
//                     if (yInterleaved)
//                     {
//                         cout << " smally: " << smally;
//                     }
//                     cout << endl;

                    if (bigxStart == bigxThresh)
                        xInterleaved = true;
                    else
                        xInterleaved = false;

                    smallx = initsmallx;                    
                    
                    PidxBY = (z * ydim + bigy) * 2;
                    if (P[PidxBY] != -1)
                    {
                        DptrBY = D + P[PidxBY];
                        if (bigxInc == 1)
                        {
                            DptrBYBX = DptrBY;
                            DptrBYSX = DptrBY + P[PidxBY + 1] - 1;
                        }
                        else
                        {
                            DptrBYBX = DptrBY + P[PidxBY + 1] - 1;
                            DptrBYSX = DptrBY;
                        }
                    }
                    else
                    {
                        DptrBY = (ShellVoxel*)NULL;
                    }

                    // FIXME: only necessary if yInterleaved (see yin == 0)
                    PidxSY = (z * ydim + smally) * 2;
                    if (P[PidxSY] != -1)
                    {
                        DptrSY = D + P[PidxSY];
                        if (bigxInc == 1)
                        {
                            DptrSYBX = DptrSY;
                            DptrSYSX = DptrSY + P[PidxSY + 1] - 1;
                        }
                        else
                        {
                            DptrSYBX = DptrSY + P[PidxSY + 1] - 1;
                            DptrSYSX = DptrSY;
                        }
                        
                    }
                    else
                    {
                        DptrSY = (ShellVoxel*)NULL;
                    }
                    

                    for (int bigx = bigxStart; bigx != bigxEnd;
                         bigx += bigxInc)
                    {

                        // -----------
                        
                        if (DptrBY)
                        {
                            if (DptrBYBX->x == bigx)
                            {
                                // render BYBX here
                                DrawVoxel(DptrBYBX, octantIdxBYBX, bigy, z,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrBYBX += bigxInc;
                            }
                            
                            if (xInterleaved && DptrBYSX->x == smallx)
                            {
                                // render BYSX here
                                DrawVoxel(DptrBYSX, octantIdxBYSX, bigy, z,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrBYSX -= bigxInc;
                            }
                        }

                        if (yInterleaved && DptrSY)
                        {
                            if (DptrSYBX->x == bigx)
                            {
                                // render SYBX here
                                DrawVoxel(DptrSYBX, octantIdxSYBX, smally, z,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrSYBX += bigxInc;
                            }

                            if (xInterleaved && DptrSYSX->x == smallx)
                            {
                                DrawVoxel(DptrSYSX, octantIdxSYSX, smally, z,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrSYSX -= bigxInc;
                            }
                        }


                        // -----------


                        if (xInterleaved)
                        {
                            smallx -= bigxInc;
                        }
                        else
                        {
                            if (bigx + bigxInc == bigxThresh)
                                xInterleaved = true;
                        }
                        
                    } // for (int bigx ...

                    if (yInterleaved)
                    {
                        smally -= bigyInc;
                    }
                    else
                    {
                        if (bigy + bigyInc == bigyThresh)
                            yInterleaved = true;
                    }
                    
                } // for (int bigy ...
            } // for (int z ...
            cout << "END: smally == " << smally << endl;            
         } // if (zin == 0) ...


         
         else if (yin == 0)
         {
             // our viewpoint is above or below the volume (in our conventional
             // coordinate system)
             
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;

            cout << "HINK: face-on, yin == 0, pq == "
                 << pq << " ps == " << ps << endl;


            // y is th slowest changing dimension
            int y0, y1, yinc;
            unsigned char octantIdxBZBX, octantIdxBZSX;
            unsigned char octantIdxSZBX, octantIdxSZSX;
            
            
            if (camVoxelPos[1] >= ydim)
            {
                y0 = 0;
                y1 = ydim;
                yinc = 1;
                // zyx == 010b == 0x2
                octantIdxBZBX =
                    octantIdxBZSX = octantIdxSZBX = octantIdxSZSX = 0x2;
            }
            else
            {
                y0 = ydim - 1;
                y1 = -1; // we're going to use != y1 as end cond
                yinc = -1;
                octantIdxBZBX =
                    octantIdxBZSX = octantIdxSZBX = octantIdxSZSX = 0;
            }

            // setup split-loop for z
            int bigzStart, bigzEnd, bigzInc, bigzThresh, smallz, initsmallz;
            if (ps >= (zdim - ps))
            {
                bigzStart = 0;
                bigzEnd = ps;
                bigzInc = 1;
                bigzThresh = ps - (zdim - ps); // biglen - small_len
                initsmallz = zdim - 1;
                octantIdxBZBX |= 0x04; octantIdxBZSX = octantIdxBZBX;
            }
            else
            {
                bigzStart = zdim - 1;
                bigzEnd = ps - 1; // loop uses != bigzEnd
                bigzInc = -1;
                bigzThresh = ps * 2 - 1;
                initsmallz = 0;
                octantIdxSZBX |= 0x4; octantIdxSZSX = octantIdxSZBX;
            }
            

            // setup split-loop for x
            int bigxStart, bigxEnd, bigxInc, bigxThresh, smallx, initsmallx;
            if (pq >= (xdim - pq))
            {
                bigxStart = 0;
                bigxEnd = pq;
                bigxInc = 1;
                bigxThresh = pq - (xdim - pq); // biglen - small_len
                initsmallx = xdim - 1;
                octantIdxBZBX |= 0x1; octantIdxSZBX = octantIdxBZBX;
            }
            else
            {
                bigxStart = xdim - 1;
                bigxEnd = pq - 1; // loop uses != bigxEnd
                bigxInc = -1;
                bigxThresh = pq * 2 - 1;
                initsmallx = 0;
                octantIdxBZSX |= 0x1; octantIdxSZSX = octantIdxBZSX;
            }

            bool zInterleaved;
            bool xInterleaved;

            ShellVoxel *DptrBZ, *DptrBZBX, *DptrBZSX;
            DptrBZ = DptrBZBX = DptrBZSX = (ShellVoxel*)NULL;
            
            ShellVoxel *DptrSZ, *DptrSZBX, *DptrSZSX;
            DptrSZ = DptrSZBX = DptrSZSX = (ShellVoxel*)NULL;
            
            int PidxBZ, PidxSZ;

            for (int y = y0; y != y1; y += yinc)
            {
                if (bigzStart == bigzThresh)
                    zInterleaved = true;
                else
                    zInterleaved = false;

                smallz = initsmallz;

                for (int bigz = bigzStart; bigz != bigzEnd; bigz += bigzInc)
                {
                    if (bigxStart == bigxThresh)
                        xInterleaved = true;
                    else
                        xInterleaved = false;

                    smallx = initsmallx;

                    PidxBZ = (bigz * ydim + y) * 2;
                    if (P[PidxBZ] != -1)
                    {
                        // ?? if !zInterleaved, we don't get here either!
                        DptrBZ = D + P[PidxBZ];
                        if (bigxInc == 1)
                        {
                            // bigx is increasing, this means that the Dptr
                            // for bigx is at the start of the sparse X
                            // structure
                            DptrBZBX = DptrBZ;
                            DptrBZSX = DptrBZ + P[PidxBZ + 1] - 1;
                        }
                        else
                        {
                            DptrBZBX = DptrBZ + P[PidxBZ + 1] - 1;
                            DptrBZSX = DptrBZ;
                        }
                    }
                    else
                    {
                        DptrBZ = (ShellVoxel*)NULL;
                    }

                    if (zInterleaved)
                    {
                        PidxSZ = (smallz * ydim + y) * 2;
                        if (P[PidxSZ] != -1)
                        {
                            DptrSZ = D + P[PidxSZ];
                            if (bigxInc == 1)
                            {
                                DptrSZBX = DptrSZ;
                                DptrSZSX = DptrSZ + P[PidxSZ + 1] - 1;
                            }
                            else
                            {
                                DptrSZBX = DptrSZ + P[PidxSZ + 1] - 1;
                                DptrSZSX = DptrSZ;
                            }
                        }
                        else
                        {
                            DptrSZ = (ShellVoxel*)NULL;
                        }
                    }
                    else
                    {
                        // zInterleaved == false, so we set this for safety
                        DptrSZ = (ShellVoxel*)NULL;
                    }

                    for (int bigx = bigxStart; bigx != bigxEnd;
                         bigx += bigxInc)
                    {

                        if (DptrBZ)
                        {
                            // ?? if !zInterleaved, we never get here!!
                            if (DptrBZBX->x == bigx)
                            {
                                DrawVoxel(DptrBZBX, octantIdxBZBX, y, bigz,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrBZBX += bigxInc;
                            }

                            if (xInterleaved && DptrBZSX->x == smallx)
                            {
                                DrawVoxel(DptrBZSX, octantIdxBZSX, y, bigz,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrBZSX -= bigxInc;
                            }
                        } // if (DptrBZ) ...

                        if (zInterleaved && DptrSZ)
                        {
                            if (DptrSZBX->x == bigx)
                            {
                                DrawVoxel(DptrSZBX, octantIdxSZBX, y, smallz,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrSZBX += bigxInc;
                            }

                            if (xInterleaved && DptrSZSX->x == smallx)
                            {
                                DrawVoxel(DptrSZSX, octantIdxSZSX, y, smallz,
                                          prev_colour, ambient, diffuse, u, v);
                                DptrSZSX -= bigxInc;
                            }
                        } // if (zInterleaved && DptrSZ) ...


                        if (xInterleaved)
                        {
                            smallx -= bigxInc;
                        }
                        else
                        {
                            if (bigx + bigxInc == bigxThresh)
                                xInterleaved = true;
                        }
                        
                        
                    } // for (int bigx == bigxStart ...

                    if (zInterleaved)
                    {
                        smallz -= bigzInc;
                    }
                    else
                    {
                        if (bigz + bigzInc == bigzThresh)
                            zInterleaved = true;
                    }
                    
                } // for (int bigz == bigzStart ...
                
            } // for (int y = y0 ...
            
         } // else if (yin == 0) ...
         
         else // xin == 0
         {
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;

            cout << "HINK: face-on, xin == 0, pp == "
                 << pr << " ps == " << ps << endl;


            // x is the slowest changing dimension
            int x0, x1, xinc;
            unsigned char octantIdxBZBY, octantIdxBZSY;
            unsigned char octantIdxSZBY, octantIdxSZSY;
            
            
            if (camVoxelPos[0] >= xdim)
            {
                x0 = 0;
                x1 = xdim;
                xinc = 1;
                // zyx == 001b == 0x1
                octantIdxBZBY = 
                    octantIdxBZSY = octantIdxSZBY = octantIdxSZSY = 0x1;
            }
            else
            {
                x0 = xdim - 1;
                x1 = -1; // we're going to use != x1 as end cond
                xinc = -1;
                octantIdxBZBY =
                    octantIdxBZSY = octantIdxSZBY = octantIdxSZSY = 0;
            }


            // setup split-loop for z
            int bigzStart, bigzEnd, bigzInc, bigzThresh, smallz, initsmallz;
            if (ps >= (zdim - ps))
            {
                bigzStart = 0;
                bigzEnd = ps;
                bigzInc = 1;
                bigzThresh = ps - (zdim - ps); // biglen - small_len
                initsmallz = zdim - 1;
                octantIdxBZBY |= 0x04; octantIdxBZSY = octantIdxBZBY;
            }
            else
            {
                bigzStart = zdim - 1;
                bigzEnd = ps - 1; // loop uses != bigzEnd
                bigzInc = -1;
                bigzThresh = ps * 2 - 1;
                initsmallz = 0;
                octantIdxSZBY |= 0x4; octantIdxSZSY = octantIdxSZBY;
            }
            
            // setup split-loop for y
            int bigyStart, bigyEnd, bigyInc, bigyThresh, smally, initsmally;
            if (pr >= (ydim - pr))
            {
                bigyStart = 0;
                bigyEnd = pr;
                bigyInc = 1;
                bigyThresh = pr - (ydim - pr); // biglen - small_len
                initsmally = ydim - 1;
                octantIdxBZBY |= 0x2; octantIdxSZBY = octantIdxBZBY;
            }
            else
            {
                bigyStart = ydim - 1;
                bigyEnd = pr - 1; // loop uses != bigyEnd
                bigyInc = -1;
                bigyThresh = pr * 2 - 1;
                initsmally = 0;
                octantIdxBZSY |= 0x2; octantIdxSZSY = octantIdxBZSY;
            }

            bool zInterleaved;
            bool yInterleaved;

            // now setup large Dptr matrix... x iteratates slowly,
            // for each new x, we might have to render something from
            // a max of 4 Dptr.  Each time we render, that specific Dptr
            // is incremented/decremented depending on the direction of x
            ShellVoxel **DptrMatrix;
            DptrMatrix = new ShellVoxel*[zdim * ydim];
            int *xMatrix = new int[zdim * ydim];

            int curIdx, curPidx;
            for (int z = 0; z < zdim; z++)
            {
                for (int y = 0; y < ydim; y++)
                {
                    curIdx = z * ydim + y;
                    curPidx = 2 * curIdx;
                    
                    if (P[curPidx] != -1)
                    {
                        DptrMatrix[curIdx] = D + P[curPidx];
                        if (xinc == -1)
                        {
                            // if x is moving "backwards", we have to
                            // initialize to the Dptr at the end
                            DptrMatrix[curIdx] += (P[curPidx + 1] - 1);
                        }
                        xMatrix[curIdx] = DptrMatrix[curIdx]->x;
                    }
                    else
                    {
                        DptrMatrix[curIdx] = (ShellVoxel*)NULL;
                        xMatrix[curIdx] = -1;
                    }
                } // for (int y = 0...
            } // for (int z = 0...

            int curBZBYidx, curBZSYidx, curSZBYidx, curSZSYidx;
            ShellVoxel **curDpp;
            int *curxPtr;

            for (int x = x0; x != x1; x += xinc)
            {
                if (bigzStart == bigzThresh)
                    zInterleaved = true;
                else
                    zInterleaved = false;

                smallz = initsmallz;

                for (int bigz = bigzStart; bigz != bigzEnd; bigz += bigzInc)
                {
                    if (bigyStart == bigyThresh)
                        yInterleaved = true;
                    else
                        yInterleaved = false;

                    smally = initsmally;

                    // we can start setting this up here (saves time in the
                    // inner loop)
                    //curBZBYidx = curBZSYidx = bigz * ydim;
                    //curSZBYidx = curSZSYidx = smallz * ydim;
                    
                    for (int bigy = bigyStart; bigy != bigyEnd;
                         bigy += bigyInc)
                    {

                        curBZBYidx = bigz * ydim + bigy;
                        curDpp = DptrMatrix + curBZBYidx;
                        curxPtr = xMatrix + curBZBYidx;

                        // we want to avoid fetching (**curDpp) EVERY time,
                        // i.e. *curDpp)->x == x would be DEADLY (fetch whole
                        // ShellVoxel struct for EVERY ITERATION
                        // *curxPtr will be -1 if (*curDpp) == NULL, so that's
                        // one check less
                        if (*curxPtr == x)
                        {
                            DrawVoxel(*curDpp, octantIdxBZBY, bigy, bigz,
                                      prev_colour, ambient, diffuse, u, v);
                            // advance to next ShellVoxel*
                            *curDpp += xinc;
                            // also perform caching of x (this will cause
                            // another lookup, but seeing that we've just
                            // done this anyway, the prefetch might have
                            // gotten it
                            *curxPtr = (*curDpp)->x;
                            //DptrMatrix[curBZBYidx] += xinc;
                            //DptrMatrix[curBZBYidx] = curDptr + xinc;
                        }

                        if (zInterleaved)
                        {
                            curSZBYidx = smallz * ydim + bigy;
                            curDpp = DptrMatrix + curSZBYidx;
                            curxPtr = xMatrix + curSZBYidx;

                            if (*curxPtr == x)
                            {
                                DrawVoxel(*curDpp, octantIdxSZBY, bigy, smallz,
                                          prev_colour, ambient, diffuse, u, v);
                                *curDpp += xinc;
                                *curxPtr = (*curDpp)->x;
                            }

                            
                            if (yInterleaved)
                            {
                                curSZSYidx = smallz * ydim + smally;
                                curDpp = DptrMatrix + curSZSYidx;
                                curxPtr = xMatrix + curSZSYidx;

                                if (*curxPtr == x)
                                {
                                    DrawVoxel(*curDpp, octantIdxSZSY,
                                              smally, smallz, prev_colour,
                                              ambient, diffuse, u, v);
                                    *curDpp += xinc;
                                    *curxPtr = (*curDpp)->x;
                                }
                            }

                        } // if (zInterleaved)...
                        
                        if (yInterleaved)
                        {
                            curBZSYidx = bigz * ydim + smally;
                            curDpp = DptrMatrix + curBZSYidx;
                            curxPtr = xMatrix + curBZSYidx;

                            if (*curxPtr == x)
                            {
                                DrawVoxel(*curDpp, octantIdxBZSY, smally, bigz,
                                          prev_colour, ambient, diffuse, u, v);
                                *curDpp += xinc;
                                *curxPtr = (*curDpp)->x;
                            }

                            smally -= bigyInc;
                        }
                        else
                        {
                            // check if next y-loop will be interleaved
                            if (bigy + bigyInc == bigyThresh)
                                yInterleaved = true;
                        }

                        
                    } // for (int bigy == bigyStart ... INNER LOOP

                    if (zInterleaved)
                    {
                        smallz -= bigzInc;
                    }
                    else
                    {
                        if (bigz + bigzInc == bigzThresh)
                            zInterleaved = true;
                    }
                    
                } // for (int bigz == bigzStart ...
                
            } // for (int x = x0 ...

            // get rid of the DptrMatrix
            delete DptrMatrix;
            delete xMatrix;
            
         } // else (xin == 0 face-on case)
      } // FACE-ON

      // edge-on
      else if (xyztot == 1)
      {

          // for edge-on, we know that the slowest changing dimension should
          // be any of the non-edge (i.e. NOT the IN) dimensions, e.g. in
          // the case of zin, x or y should preferably be the slowest changing.
          // because the two fastest changing ones are always angled away
          // from us, their order with respect to each other is not critical
          
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "EDGE-ON" << endl;
#endif         
         int edge_idx = 0;
         if (zin)
         {

            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;
            cout << "HONK: edge-on, zin, ps == " << ps << endl;            

            int x0, x1, xinc, y0, y1, yinc;            
            if (camVoxelPos[0] >= xdim)
            {
                x0 = 0;
                x1 = xdim;
                xinc = 1;
                edge_idx |= 0x1;                
            }
            else
            {
                x0 = xdim - 1;
                x1 = -1; // we're going to use != x1 as end condition
                xinc = -1;
            }
            if (camVoxelPos[1] >= ydim)
            {
                y0 = 0;
                y1 = ydim;
                yinc = 1;
               edge_idx |= 0x2;
            }
            else
            {
                y0 = ydim - 1;
                y1 = -1;
                yinc = -1;
            }

	    // start of experimental iPBTF code -----------------------------
	    // hmmm... split outer loop
	    int seg1len = ps;
            int seg2len = zdim - ps;

            int bigzStart, bigzEnd, bigzInc, bigzThresh, smallz, initsmallz;
            unsigned char bigOctantIdx, smallOctantIdx;
            if (seg1len >= seg2len)
            {
                initsmallz = zdim - 1;
                bigzStart = 0;
                bigzEnd = ps;
                bigzInc = 1;
                // if bigz >= bigThresh, it means we have to interleave
                bigzThresh = seg1len - seg2len;
                bigOctantIdx = edge_idx + 4;
                smallOctantIdx = edge_idx;
            }
            else
            {
                initsmallz = 0;
                bigzStart = zdim - 1;
                bigzEnd = ps - 1; // our loop is going to use !=
                bigzInc = -1;
                bigzThresh = seg1len * 2 - 1;
                bigOctantIdx = edge_idx;
                smallOctantIdx = edge_idx + 4;
            }

            
            bool interleaving;

            

            for (int y = y0; y != y1; y += yinc)
            {
                //for (int bigz = bigzStart; bigz != bigzEnd; bigz += bigzInc)
                //{
//                 if (bigz == smallz)
//                 {
//                     cerr << "AAAARRRGGGHHH THIS SHOULD NEVER HAPPEN!" << endl;
//                 }

                smallz = initsmallz;
                if (bigzStart == bigzThresh)
                    interleaving = true;
                else
                    interleaving = false;
                
                for (int bigz = bigzStart; bigz != bigzEnd;
                     bigz += bigzInc)
                {                

                    if (interleaving)
                    {
                        Pidx0 = (bigz * ydim + y) * 2;
                        Pidx1 = (smallz * ydim + y) * 2;
                        
                        int len0 = 0, len1 = 0;
                        
                        if (P[Pidx0] != -1)
                        {
                            Dptr0 = D + P[Pidx0];
                            len0 = P[Pidx0 + 1];
                            if (xinc == -1)
                                Dptr0 += (len0 - 1);
                        }
                        if (P[Pidx1] != -1)
                        {
                            Dptr1 = D + P[Pidx1];
                            len1 = P[Pidx1 + 1];
                            if (xinc == -1)
                                Dptr1 += (len1 - 1);
                        }

                        if (len0 && len0 >= len1)
                        {
                            for (int x = 0; x < len0; x++)
                            {
                                DrawVoxel(Dptr0, bigOctantIdx, y, bigz,
                                          prev_colour, ambient, diffuse, u, v);
                                Dptr0 += xinc;
                                if (x < len1)
                                {
                                    DrawVoxel(Dptr1, smallOctantIdx, y, smallz,
                                              prev_colour, ambient, diffuse,
                                              u, v);
                                    Dptr1 += xinc;
                                }
                            }
                        }

                        if (len1 && len1 > len0)
                        {
                            for (int x = 0; x < len1; x++)
                            {
                                DrawVoxel(Dptr1, smallOctantIdx, y, smallz,
                                          prev_colour, ambient, diffuse, u, v);
                                Dptr1 += xinc;
                                
                                if (x < len0)
                                {
                                    DrawVoxel(Dptr0, bigOctantIdx, y, bigz,
                                              prev_colour, ambient, diffuse,
                                              u, v);
                                    Dptr0 += xinc;
                                }
                            }
                        }

                        // we're interleaving, so increment smallz (opposite dir)
                        smallz -= bigzInc;
                        
                        
                    } // if (interleaving) ...

                    else
                    {
                        Pidx = (bigz * ydim + y) * 2;
                        if (P[Pidx] != -1)
                        {
                            Dptr = D + P[Pidx];
                            if (xinc == -1)
                                Dptr += (P[Pidx + 1] - 1);

                            for (int x = 0; x < P[Pidx + 1]; x++)
                            {
                                DrawVoxel(Dptr, bigOctantIdx, y, bigz,
                                          prev_colour, ambient, diffuse, u, v);
                                Dptr += xinc;
                            } // for (int x ...
                        } // if (P[Pidx] != -1 ...
                    

                        // see whether we have to interleave on the next loop
                        if (bigz + bigzInc == bigzThresh)
                        {
                            interleaving = true;
                        }
                    } // else interleaving ...
                } // for (int bigz ...
            } // for (int y = ...
	    
         } // if (zin ...
         else if (yin)
         {
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            if (camVoxelPos[0] >= xdim)
               edge_idx |= 0x1;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, xdim, 0, pr + 1, 0, zdim, edge_idx + 2, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(0, xdim, pr + 1, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
         else // xin
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            if (camVoxelPos[1] >= ydim)
               edge_idx |= 0x2;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, pq + 1,      0, ydim, 0, zdim, edge_idx + 1, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(pq + 1, xdim,   0, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
      } // EDGE-ON

      // corner-on
      else
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "CORNER-ON" << endl;
#endif         
         // being here means Vp is by definition outside the volume:
         octantIdx = 0;
         if (camVoxelPos[0] >= inputDims[0])
         {
            octantIdx |= 0x1;
         }
         if (camVoxelPos[1] >= inputDims[1])
         {
            octantIdx |= 0x2;
         }
         if (camVoxelPos[2] >= inputDims[2])
         {
            octantIdx |= 0x4;
         }
         this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2], octantIdx, D, P, inputDims[1], ambient, diffuse, u, v);
      } // CORNER-ON

   } // else PERSPECTIVE rendering with interleaved PBTF

   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   else  // PERSPECTIVE rendering, extra-special BTF (mode 2)
   {

       // we need to project the centre of the volume onto the view plane and
       // then check in which volume octant that point is

       vtkMatrix4x4 *octantM = vtkMatrix4x4::New();
       vtkMatrix4x4::Multiply4x4(this->PerspectiveMatrix,
                                 this->VoxelsToViewMatrix,  octantM);
       float voxel_volume_centre[4];
       float projected_volume_centre[4];
       float projected_voxel_volume_centre[4];
       
       voxel_volume_centre[0] = (float)xdim / 2.0;
       voxel_volume_centre[1] = (float)ydim / 2.0;
       voxel_volume_centre[2] = (float)zdim / 2.0;
       voxel_volume_centre[3] = 1.0;

       // project the volume centre onto the view plane
       octantM->MultiplyPoint(voxel_volume_centre, projected_volume_centre);

       // set its third dimension to zero (that should be IN the plane)
       projected_volume_centre[2] = 0.0;

       // we're going back to voxel space
       octantM->Invert();
       octantM->MultiplyPoint(projected_volume_centre,
                              projected_voxel_volume_centre);

       // now determine the octant
       // and yes, it does make sense that we're checking with <: the projected
       // point is not stuck on the view-plane, it's slipping and sliding as
       // the view plane is moving.  Think on that and take it from there. :)
       octantIdx = 0;
       if (projected_voxel_volume_centre[0] < voxel_volume_centre[0])
       {
           octantIdx |= 0x1;
       }
       if (projected_voxel_volume_centre[1] < voxel_volume_centre[1])
       {
           octantIdx |= 0x2;
       }
       if (projected_voxel_volume_centre[2] < voxel_volume_centre[2])
       {
           octantIdx |= 0x4;
       }

      // take care of the matrix we made
      octantM->Delete();

#ifdef SSM_VERBOSE_OUTPUT
      cout << "PVC " << projected_voxel_volume_centre[0] << " " << projected_voxel_volume_centre[1] << " " << projected_voxel_volume_centre[2] << " " << endl;
#endif

      // check whether it's volume-on, face-on, edge-on or corner-on
      int xin = 0, yin = 0, zin = 0, xyztot = 0;
      if (camVoxelPos[0] >= 0 && camVoxelPos[0] < inputDims[0])
      {
         xin = 1;
         xyztot++;
      }
      if (camVoxelPos[1] >= 0 && camVoxelPos[1] < inputDims[1])
      {
         yin = 1;
         xyztot++;
      }
      if (camVoxelPos[2] >= 0 && camVoxelPos[2] < inputDims[2])
      {
         zin = 1;
         xyztot++;
      }

      int pq, pr, ps;

      // **************************************************************************
      // NB: convention: LDB is our origin, NOT LDF as in Swan's thesis!
      // **************************************************************************

      // volume-on: break up volume into 8 octants
      if (xyztot == 3)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "VOLUME-ON" << endl;
#endif         
         pq = vtkMath::Round(camVoxelPos[0]);
         pq = (pq >= xdim) ? xdim - 1 : pq;
         pr = vtkMath::Round(camVoxelPos[1]);
         pr = (pr >= ydim) ? ydim - 1 : pr;
         ps = vtkMath::Round(camVoxelPos[2]);
         ps = (ps >= zdim) ? zdim - 1 : ps;
         // rendering octant LDF (code: 4)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      ps + 1, zdim,   3, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDF (code: 5)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      ps + 1, zdim,   2, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUF (code: 6)      
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   ps + 1, zdim,   1, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUF (code: 7)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   ps + 1, zdim,   0, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LDB (code: 0)
         this->DrawVoxels(0, pq + 1,      0, pr + 1,      0, ps + 1,      7, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RDB (code: 1)
         this->DrawVoxels(pq + 1, xdim,   0, pr + 1,      0, ps + 1,      6, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant LUB (code: 2)
         this->DrawVoxels(0, pq + 1,      pr + 1, ydim,   0, ps + 1,      5, D, P, ydim, ambient, diffuse, u, v);
         // rendering octant RUB (code: 3)
         this->DrawVoxels(pq + 1, xdim,   pr + 1, ydim,   0, ps + 1,      4, D, P, ydim, ambient, diffuse, u, v);
      } // end VOLUME-ON

      // face-on: break volume into 4 quadrants
      else if (xyztot == 2)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "FACE-ON: ";
#endif
         // check in which octant pq,pr,ps falls
         octantIdx = 0;

         // this is face-on, only one dim can be outside
         if (zin == 0)
         {
             // z is outside
             if (camVoxelPos[2] >= zdim)
             {
                 octantIdx |= 0x4;
             }
             if (camVoxelPos[0] > xdim - camVoxelPos[0])
                 octantIdx |= 0x1;
             if (camVoxelPos[1] > ydim - camVoxelPos[1])
                 octantIdx |= 0x2;
         }
         else if (yin == 0)
         {
             // y is outside
             if (camVoxelPos[1] >= ydim)
             {
                 octantIdx |= 0x2;
             }
             if (camVoxelPos[0] > xdim - camVoxelPos[0])
                 octantIdx |= 0x1;
             if (camVoxelPos[2] > zdim - camVoxelPos[2])
                 octantIdx |= 0x4;
             
         }
         else if (xin == 0)
         {
             // x is outside
             if (camVoxelPos[0] >= xdim)
             {
                 octantIdx |= 0x1;
             }
             if (camVoxelPos[1] > ydim - camVoxelPos[1])
                 octantIdx |= 0x2;
             if (camVoxelPos[2] > zdim - camVoxelPos[2])
                 octantIdx |= 0x4;
         }
             
         this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2],
                          octantIdx, D, P, inputDims[1], ambient, diffuse,
                          u, v);
         
      } // FACE-ON

      // edge-on
      else if (xyztot == 1)
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "EDGE-ON" << endl;
#endif         
         int edge_idx = 0;
         if (zin)
         {
            ps = vtkMath::Round(camVoxelPos[2]);
            ps = (ps >= zdim) ? zdim - 1 : ps;
            if (camVoxelPos[0] >= xdim)
               edge_idx |= 0x1;
            if (camVoxelPos[1] >= ydim)
               edge_idx |= 0x2;
            // this is sorted, really
            this->DrawVoxels(0, xdim,   0, ydim,   0, ps + 1,      edge_idx + 4, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(0, xdim,   0, ydim,   ps + 1, zdim,   edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
         else if (yin)
         {
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            if (camVoxelPos[0] >= xdim)
               edge_idx |= 0x1;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, xdim, 0, pr + 1, 0, zdim, edge_idx + 2, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(0, xdim, pr + 1, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
         else // xin
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            if (camVoxelPos[1] >= ydim)
               edge_idx |= 0x2;
            if (camVoxelPos[2] >= zdim)
               edge_idx |= 0x4;
            this->DrawVoxels(0, pq + 1,      0, ydim, 0, zdim, edge_idx + 1, D, P, inputDims[1], ambient, diffuse, u, v);
            this->DrawVoxels(pq + 1, xdim,   0, ydim, 0, zdim, edge_idx, D, P, inputDims[1], ambient, diffuse, u, v);
         }
      } // EDGE-ON

      // corner-on
      else
      {
#ifdef SSM_VERBOSE_OUTPUT          
         cout << "CORNER-ON" << endl;
#endif         
         // being here means Vp is by definition outside the volume:
         octantIdx = 0;
         if (camVoxelPos[0] >= inputDims[0])
         {
            octantIdx |= 0x1;
         }
         if (camVoxelPos[1] >= inputDims[1])
         {
            octantIdx |= 0x2;
         }
         if (camVoxelPos[2] >= inputDims[2])
         {
            octantIdx |= 0x4;
         }
         this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2], octantIdx, D, P, inputDims[1], ambient, diffuse, u, v);
      } // CORNER-ON

   } // else PERSPECTIVE rendering, extra-special BTF
   
#ifdef SSM_VERBOSE_OUTPUT
   cout << "voxels_drawn == " << voxels_drawn << endl;
#endif   

   // end plotting frikking points
   glEnd();

   if (RenderMode == 0)
   {
     glDisable(GL_TEXTURE_2D);
#ifdef GL_VERSION_1_1
     GLuint tempIndex = tex_idx;
     glDeleteTextures(1, &tempIndex);
#else
     glDeleteLists(tex_idx, 1);
#endif
   }

   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();

   // we have to reactivate this, or the VTK z-clear doesn't work!
   glDepthMask(GL_TRUE);

   // let's deactivate blending again
   glDisable(GL_BLEND);

#ifdef SSM_VERBOSE_OUTPUT
   clock_t end_clock = clock();
   clock_t diff_clock = end_clock - start_clock;
   
   float secs = (float)diff_clock / (float)CLOCKS_PER_SEC;   
   cout << "Clock ticks == " << diff_clock << " Secs == " << secs << " FPS == " << 1.0 / secs << endl;
#endif   

   // restore all gl attributes
   glPopAttrib();
}


// we fill out a 3 bit octant field, the value is the octant index
// lsb is x, msb is z
// SO: 0 == 000 == x0y0z0
//       1 == 001 == x1y0z0
//       2 == 010 == x0y1z0
// etc.

void vtkOpenGLVolumeShellSplatMapper::DrawVoxels(int x0, int x1, 
						 int y0, int y1, 
                                                 int z0, int z1, 
						 unsigned char octantIdx, 
						 ShellVoxel *D, int *P, 
						 int ydim,
                                                 const float& ambient, 
						 const float& diffuse, 
						 float* u, float* v)
{

   int x, y, z;
   int Pidx;
   ShellVoxel* Dptr;
   GLfloat prev_colour[4];
   prev_colour[0] = prev_colour[1] = prev_colour[2] = prev_colour[3] = -1;
#ifdef SSM_VERBOSE_OUTPUT   
   cout << x0 << " --> " << x1 << ", " << y0 << " --> " << y1 << ", " << z0 << " --> " << z1 << endl;
#endif   
   switch (octantIdx)
   {
      case 0:
         {
            // octant 0, BTF: z-, y-, x-
            for (z = z1 - 1; z >= z0; z--)
            {
               for (y = y1 - 1; y >= y0; y--)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx] + P[Pidx + 1] - 1;
                     // if Dptr->x < x0 we can cop out here already! - should we do this?
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z,
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr--;
                        } // if (Dptr->x ...
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 1:
         {
            // octant 1, BTF: z-, y-, x+
            for (z = z1 - 1; z >= z0; z--)
            {
               for (y = y1 - 1; y >= y0; y--)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx];
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr++;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 2:
         {
            // octant 2, BTF: z-, y+, x-
            for (z = z1 - 1; z >= z0; z--)
            {
               for (y = y0; y < y1; y++)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx] + P[Pidx + 1] - 1;                        
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr--;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 3:
         {
            // octant 3, BTF: z-, y+, x+
            for (z = z1 - 1; z >= z0; z--)
            {
               for (y = y0; y < y1; y++)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx];
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr++;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 4:
         {
            // octant 4, BTF: z+, y-, x-
            for (z = z0; z < z1; z++)
            {
               for (y = y1 - 1; y >= y0; y--)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx] + P[Pidx + 1] - 1;
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr--;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 5:
         {
            // octant 5, BTF: z+, y-, x+
            for (z = z0; z < z1; z++)
            {
               for (y = y1 - 1; y >= y0; y--)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx];
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr++;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
         break;

      case 6:
         {
            // octant 6, BTF: z+, y+, x-
            for (z = z0; z < z1; z++)
            {
               for (y = y0; y < y1; y++)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx] + P[Pidx + 1] - 1;                        
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr--;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
#ifdef SSM_VERBOSE_OUTPUT         
         cout << endl;
#endif         
         break;

      case 7:
         {
            // octant 7, BTF: z+, y+, x+
            for (z = z0; z < z1; z++)
            {
               for (y = y0; y < y1; y++)
               {
                  Pidx = (z * ydim + y) * 2;
                  // only do something if we have shellvoxels for this x-row                    
                  if (P[Pidx] != -1)
                  {
                     Dptr = D + P[Pidx];                        
                     for (x = 0; x < P[Pidx + 1]; x++)
                     {
                        if (Dptr->x >= x0 && Dptr->x < x1)
                        {
                           DrawVoxel(Dptr, octantIdx, y, z, 
                                     prev_colour, ambient, diffuse, u, v);
                           Dptr++;
                        }
                     } // for (x = 0 ...
                  }
               } // for (z = inputDims[2] ...
            } // for (y = inputDims[1] ...
         }
#ifdef SSM_VERBOSE_OUTPUT         
         cout << endl;
#endif         
         break;
   }

}
