// vtkShellExtractor.h copyright (c) 2003 
// by Charl P. Botha cpbotha@ieee.org 
// and the TU Delft Visualisation Group http://visualisation.tudelft.nl/
// $Id: vtkShellExtractor.h,v 1.3 2003/09/03 11:13:13 cpbotha Exp $
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

#ifndef __vtkShellExtractor_h
#define __vtkShellExtractor_h

#include "vtkColorTransferFunction.h"
#include "vtkImageData.h"
#include "vtkObject.h"
#include "vtkPiecewiseFunction.h"
#include "vtkTimeStamp.h"
#include "vtkdscasRenderingWin32Header.h"

struct ShellVoxel
{
    /**
     * x-coordinate of voxel.
     */
    int x;
    float Value;
    float Opacity;
    float Red;
    float Green;
    float Blue;
    /**
     * 3-component normal at voxel in volume-space.  Spacings in all directions
     * assumed to be 1.0, you can correct for this with a transformation
     * matrix, or by actually adjusting the normal.
     */
    float Normal[3];
    /**
     * Neighbourhood bit code.  I have assigned the bits as follows
     * (in a right-hand cartesian system):
     *  5  4  3  2  1  0
     * +z -z +y -y +x -x
     * Bits 6 and 7 are not significant.
     */
    unsigned char nbrOCode;
    /**
     * Coordinates of this voxel in volume space.
     */
    float volCoords[3];
};

/**
 * VTK class that is able to extract shells as documented in:
 * Jayaram K. Udupa and Dewey Odhner, "Shell Rendering", IEEE Computer Graphics
 * and Applications, 1993, pp 58--67.
 */
class VTK_DSCAS_RENDERING_EXPORT vtkShellExtractor : public vtkObject
{
public:
   vtkTypeMacro(vtkShellExtractor, vtkObject);
   /**
    * Setup shell extractor with defaults.  We make use of the vtk macro
    * to make a default New method.
    */
   static vtkShellExtractor* New();
   /**
    * Set vtkImageData input.  Sets the supplied input, also updates its 
    * references.
    */
   vtkSetObjectMacro(Input, vtkImageData);
   /**
    * Get the vtkImageData input.
    */
   vtkGetObjectMacro(Input, vtkImageData);
   /**
    * Set opacity transfer function that will be used in shell extraction.
    * References on the opacity transfer function will be correctly updated.
    */
   vtkSetObjectMacro(OpacityTF, vtkPiecewiseFunction);
   /**
    * Get transfer function that is currently being used.
    */
   vtkGetObjectMacro(OpacityTF, vtkPiecewiseFunction);
   vtkSetObjectMacro(ColourTF, vtkColorTransferFunction);
   vtkGetObjectMacro(ColourTF, vtkColorTransferFunction);
   /**
    * Set optional imagedata member that will be used only for calculating
    * rendering gradient.  No checking is done on the validity of this data,
    * so it's probably your responsibility.
    */
   vtkSetObjectMacro(GradientImageData, vtkImageData);
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
   
   ShellVoxel* GetD(void) {return this->D;}
   int* GetP(void) {return this->P;}
   //vtkGetMacro(D, ShellVoxel*);
   //vtkGetMacro(P, int*);
   
   /**
    * This function coordinates all the work.
    */
   void Update(void);
   
protected:
   vtkShellExtractor();
   ~vtkShellExtractor();

   vtkImageData* Input;
   vtkPiecewiseFunction* OpacityTF;
   vtkColorTransferFunction* ColourTF;
   float OmegaL;
   float OmegaH;
   vtkImageData* GradientImageData;

   /**
    * The time at which the last shell extraction was done.
    */
   vtkTimeStamp ExtractTime;
   
   /**
    * The D data structure which is an array of shell voxels.  Internally,
    * we start this up as a vector, but the VTK wrapping doesn't understand
    * this.
    */
   ShellVoxel *D;
   
   /**
    * Y by Z array with indices into D and row lengths for quick access.
    * NOTE: for each Y,Z pos there are TWO (2) ints!
    */
   int* P;
public:
   /**
    * Lookup table for visibility with neighboor opacity codes.
    */
   static const int shell_noc_visibility_lut[8][64];
};

#endif
