import vtk

reader = vtk.vtkStructuredGridReader()
reader.SetFileName("milieux.vtk")
reader.Update()
writer = vtk.vtkXMLStructuredGridWriter()
writer.SetFileName("milieux.vts")
writer.SetInputData(reader.GetOutput())
writer.SetDataModeToBinary()  # Pour ASCII : writer.SetDataModeToAscii()
writer.SetCompressorTypeToZLib()  # Compression zlib
writer.Write()