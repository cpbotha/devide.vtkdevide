// vtkInteractorStyleTrackballActorConstrained copyright (c) 2003 Charl P. Botha <cpbotha@ieee.org>
// $Id: vtkInteractorStyleTrackballActorConstrained.h,v 1.5 2004/01/15 11:00:55 cpbotha Exp $
// interactorstyle that can constrain object manipulation to planes, lines and objects

#ifndef __vtkInteractorStyleTrackballActorConstrained_h
#define __vtkInteractorStyleTrackballActorConstrained_h

#include <vtkInteractorStyleTrackballActor.h>
#include <vtkObjectFactory.h>
#include <vtkstd/vector>
#include "vtkdevideRenderingWin32Header.h"

class vtkProp3D;

class VTK_DEVIDE_RENDERING_EXPORT vtkInteractorStyleTrackballActorConstrained : public vtkInteractorStyleTrackballActor
{
public:
   static vtkInteractorStyleTrackballActorConstrained *New();
   vtkTypeRevisionMacro(vtkInteractorStyleTrackballActorConstrained, vtkInteractorStyleTrackballActor);

   virtual void OnLeftButtonDown();
   virtual void OnMiddleButtonDown();
   virtual void OnRightButtonDown();
   
   /**
    * Add a single active prop to the internal activeprops list.  Only
    * props in this list can be manipulated.  References are NOT taken
    * out.  The method will check for duplicates and quitely do nothing
    * if theProp is already in the list.
    */
   void AddActiveProp(vtkProp3D *theProp);

   /**
    * Remove theProp from the internal ActiveProps list.
    */
   void RemoveActiveProp(vtkProp3D *theProp);

   /**
    * Clear the internal list of ActiveProps.
    */
   void RemoveAllActiveProps(void);

   /**
    * Clear the current ActiveProp list and add the passed prop to it.
    * If the parameter is NULL, just clear the list.
    */
   void SetActiveProp(vtkProp3D *theProp);

   /**	
    * Return the current InteractionProp.  This can be queried after the last
	* interaction (e.g. during an EndInteraction observer) to find out if an object
	* was actually moved around.
	*/
   vtkProp3D* GetInteractionProp(void)
   {
	   return this->InteractionProp;
   }

protected:
   vtkInteractorStyleTrackballActorConstrained();
   ~vtkInteractorStyleTrackballActorConstrained();

   /**
    * Returns true if theProp is in the internal ActiveProps list.
    */
   bool IsPropActive(vtkProp3D *theProp);

//BTX
   vtkstd::vector<vtkProp3D *> ActiveProps;
//ETX

private:
   vtkInteractorStyleTrackballActorConstrained(const vtkInteractorStyleTrackballActorConstrained&);  // Not implemented.
   void operator=(const vtkInteractorStyleTrackballActor&);  // Not implemented.  
};

#endif


