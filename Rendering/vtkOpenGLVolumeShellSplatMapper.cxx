// vtkOpenGLVolumeShellSplatMapper copyright (c) 2002 by Charl P. Botha 
// http://cpbotha.net/
// $Id: vtkOpenGLVolumeShellSplatMapper.cxx,v 1.2 2003/04/29 17:10:43 cpbotha Exp $
// vtk class for volume rendering by shell splatting

// TODO:
// Use glPushAttrib, glPopAttrib to make this class less dangerous

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
   static GLfloat vt_array[12];
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

   // this will set the Gaussian parameters and initialise the textures themselves
   this->SetRenderMode(0);
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

   SetInput(NULL);
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

   clock_t start_clock = clock();

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
   cout << "nfactor == " << nfactor << endl;
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

   long tex_idx;

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

      glEndList();

      glCallList ((GLuint) tex_idx);

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
   cout << "camWorldPos == " << camWorldPos[0] << "," << camWorldPos[1] << "," << camWorldPos[2] << endl;
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

      cout << "PVC " << projected_voxel_volume_centre[0] << " " << projected_voxel_volume_centre[1] << " " << projected_voxel_volume_centre[2] << " " << endl;

      // take care of the matrix we made
      octantM->Delete();

      this->DrawVoxels(0, inputDims[0], 0, inputDims[1], 0, inputDims[2], octantIdx, D, P, inputDims[1], ambient, diffuse, u, v);
   }
   else // PERSPECTIVE rendering
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
         cout << "VOLUME-ON" << endl;
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
         cout << "FACE-ON: ";
         if (zin == 0)
         {
            pq = vtkMath::Round(camVoxelPos[0]);
            pq = (pq >= xdim) ? xdim - 1 : pq;
            pr = vtkMath::Round(camVoxelPos[1]);
            pr = (pr >= ydim) ? ydim - 1 : pr;
            if (camVoxelPos[2] >= zdim)
            {
               cout << "z >= zdim" <<endl;
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
               cout << "z < 0" <<endl;
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
         cout << "EDGE-ON" << endl;
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
         cout << "CORNER-ON" << endl;
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

   } // else PERSPECTIVE rendering

   cout << "voxels_drawn == " << voxels_drawn << endl;

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

   clock_t end_clock = clock();
   clock_t diff_clock = end_clock - start_clock;
   float secs = (float)diff_clock / (float)CLOCKS_PER_SEC;

   cout << "Clock ticks == " << diff_clock << " Secs == " << secs << " FPS == " << 1.0 / secs << endl;
}


// we fill out a 3 bit octant field, the value is the octant index
// lsb is x, msb is z
// SO: 0 == 000 == x0y0z0
//       1 == 001 == x1y0z0
//       2 == 010 == x0y1z0
// etc.

void vtkOpenGLVolumeShellSplatMapper::DrawVoxels(int x0, int x1, int y0, int y1, 
                                                 int z0, int z1, unsigned char octantIdx, ShellVoxel *D, int *P, int ydim,
                                                 const float& ambient, const float& diffuse, float* u, float* v)
{

   int x, y, z;
   int Pidx;
   ShellVoxel* Dptr;
   GLfloat prev_colour[4];
   prev_colour[0] = prev_colour[1] = prev_colour[2] = prev_colour[3] = -1;
   cout << x0 << " --> " << x1 << ", " << y0 << " --> " << y1 << ", " << z0 << " --> " << z1 << endl;
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
         cout << endl;
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
         cout << endl;
         break;
   }

}
