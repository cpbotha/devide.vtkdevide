// vtkOpenGLVolumeShellSplatMapper copyright (c) 2002 by Charl P. Botha 
// http://cpbotha.net/
// $Id: vtkOpenGLVolumeShellSplatMapper.h,v 1.3 2003/05/06 11:34:47 cpbotha Exp $
// vtk class for volume rendering by shell splatting

#ifndef __vtkOpenGLVolumeShellSplatMapper_h
# define __vtkOpenGLVolumeShellSplatMapper_h

# include "vtkObject.h"
#include "vtkImageData.h"
# include "vtkVolumeMapper.h"
# include "vtkdscasRenderingWin32Header.h"

class vtkImageData;
class vtkMatrix4x4;
class vtkRenderer;
class vtkShellExtractor;
struct ShellVoxel;
class vtkVolume;

class VTK_DSCAS_RENDERING_EXPORT vtkOpenGLVolumeShellSplatMapper : 
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

  vtkSetMacro(EllipsoidDiameter, double);
  vtkGetMacro(EllipsoidDiameter, double);
  void SetGaussianRadialExtent(double);
  double GetGaussianRadialExtent(void) {return this->gaussian_radial_extent;}
  void SetGaussianSigma(double);
  double GetGaussianSigma(void) {return this->gaussian_sigma;}

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
};

#endif
