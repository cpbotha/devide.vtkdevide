dimensions = (192,192,192)
interval = 16

import vtk

output = vtk.vtkImageData()
output.SetOrigin(0,0,0)
output.SetSpacing(0.5,0.5,1.0)
output.SetDimensions(dimensions)
output.SetNumberOfScalarComponents(1)
output.SetScalarTypeToUnsignedChar()
output.AllocateScalars()

scalars = output.GetPointData().GetScalars()
scalars.FillComponent(0, 0) # second element is fill value

c = 0
iborders = [0,0]
jborders = [0,0]
kborders = [0,0]

for idx,borders in ((0, iborders), (1, jborders), (2, kborders)):
    b = dimensions[idx] / 4
    borders[0] = b
    borders[1] = dimensions[idx] - b

for k in xrange(kborders[0],kborders[1]):
    kOdd = (k / interval) % 2
    
    for j in xrange(jborders[0], jborders[1]):
        jOdd = (j / interval) % 2
        
        for i in xrange(iborders[0], iborders[1]):
            iOdd = (i / interval) % 2
            
            if iOdd ^ jOdd ^ kOdd:
                c = 255
            else:
                c = 254

            output.SetScalarComponentFromDouble(i,j,k,0,c)

    print("%d\n" % k)
    
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName('cube.vti')
writer.SetInput(output)
writer.Write()

