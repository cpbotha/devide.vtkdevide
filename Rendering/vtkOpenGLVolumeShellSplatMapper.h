// vtkOpenGLVolumeShellSplatMapper copyright (c) 2003 
// by Charl P. Botha cpbotha@ieee.org 
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id: vtkOpenGLVolumeShellSplatMapper.h,v 1.8 2004/01/15 11:00:55 cpbotha Exp $
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


#ifndef __vtkOpenGLVolumeShellSplatMapper_h
# define __vtkOpenGLVolumeShellSplatMapper_h

# include "vtkObject.h"
#include "vtkImageData.h"
# include "vtkVolumeMapper.h"
# include "vtkdevideRenderingWin32Header.h"

class vtkImageData;
class vtkMatrix4x4;
class vtkRenderer;
class vtkShellExtractor;
struct ShellVoxel;
class vtkVolume;

class VTK_DEVIDE_RENDERING_EXPORT vtkOpenGLVolumeShellSplatMapper : 
public vtkVolumeMapper
{
 public:
  static vtkOpenGLVolumeShellSplatMapper* New();   
  vtkTypeMacro(vtkOpenGLVolumeShellSplatMapper, vtkVolumeMapper);
  /**
   * Set lower opacity bound for shell extraction.  This also calls
   * this->Modified() to indicate that the current output is invalidated.
   * Remember: 0.0 < OmegaL <= OmegaH < 1.0.
   */
  vtkSetClampMacro(OmegaL, float, 0.0, 1.0);
  /**
   * Gets current lower opacity bound for shell extraction.
   */
  vtkGetMacro(OmegaL, float);
  /**
   * Sets upper opacity bound for shell extraction.  This also calls
   * this->Modified() to indicate that the current output is invalidated.
   * Remember: 0.0 < OmegaL <= OmegaH < 1.0.
   */
  vtkSetClampMacro(OmegaH, float, 0.0, 1.0);
  /**
   * Gets current upper opacity bound for shell extraction.
   */
  vtkGetMacro(OmegaH, float);
  //vtkSetClampMacro(RenderMode, int, 0, 1);
  /**
   * Sets current RenderMode. 0 == high quality rendering with textured
   * quads, 1 == fast rendering with untextured quads, 2 == fastest rendering
   * with points.
   */
  void SetRenderMode(int newRenderMode);
  vtkGetMacro(RenderMode, int);

  /**
   * Sets current perspective ordering mode.  0 == PBTF.  1 == interleaved
   * PBTF, 2 == extra-special BTF
   */
  vtkSetClampMacro(PerspectiveOrderingMode, int, 0, 2);
  vtkGetMacro(PerspectiveOrderingMode, int);

  vtkSetMacro(EllipsoidDiameter, double);
  vtkGetMacro(EllipsoidDiameter, double);
  void SetGaussianRadialExtent(double);
  double GetGaussianRadialExtent(void) {return this->gaussian_radial_extent;}
  void SetGaussianSigma(double);
  double GetGaussianSigma(void) {return this->gaussian_sigma;}

  /**
   * If you know what yer doing, you can use this instance of the ShellExtractor
   * to do some fine tuning.
   */
  vtkGetObjectMacro(ShellExtractor, vtkShellExtractor);

  /**
   * If this is set, the gradient that's used for rendering will be calculated
   * using this vtkImageData.  Use this when you're rendering binary segmented
   * data but want to use the gradient from the original unsegmented data for
   * higher quality volume renderings.
   */
  vtkSetObjectMacro(GradientImageData, vtkImageData);

  /**
   * This gets called every time somebody wants us to draw the whole volume.
   * I'm _hoping_ this is only when viewpoint changes.
   */
  virtual void Render(vtkRenderer* ren, vtkVolume* vol);
   
 protected:
  vtkOpenGLVolumeShellSplatMapper();
  ~vtkOpenGLVolumeShellSplatMapper();

  /**
   * Internal method used by Render() to render sub-volumes.
   */
  void DrawVoxels(int x0, int x1, int y0, int y1, 
                  int z0, int z1, unsigned char octantIdx, ShellVoxel *D, int *P, int ydim,
                  const float& ambient, const float& diffuse, float* u, float* v);

  vtkMatrix4x4* volMatrix;
  vtkMatrix4x4* VoxelsToWorldMatrix;
  vtkMatrix4x4* WorldToVoxelsMatrix;
  vtkMatrix4x4* ViewMatrix;
  vtkMatrix4x4* VoxelsToViewMatrix;
  vtkMatrix4x4* ViewToVoxelsMatrix;
  /// camera to viewport matrix
  vtkMatrix4x4* PerspectiveMatrix;
   
  // 1.0 means a voxel-sized sphere will be used to create the ellipsoid
  // projections
  double EllipsoidDiameter;
  // Gaussian reconstruction function radial extent
  double gaussian_radial_extent;
  double gaussian_sigma;

  double* integrated_rfunc;
  float* normalised_integrated_rfunc;
   
  void CalculatePerViewMatrices(vtkRenderer* ren, vtkVolume* vol);   
   
  vtkShellExtractor* ShellExtractor;
  float OmegaL;
  float OmegaH;
  vtkImageData* GradientImageData;
  /**
   * Different rendering modes.  0 is for splatting, 1 is for voxel point 
   * projection.
   */
  int RenderMode;
  /**
   * Different perspective ordering modes.  0 is for Ed Swan's PBTF,
   * 1 is for interleaved PBTF (which will hopefully solve all of PBTF's
   * problems when points are not infinitesimally small)
   */
  int PerspectiveOrderingMode;
};

#endif
