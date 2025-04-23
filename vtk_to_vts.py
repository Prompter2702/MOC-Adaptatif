import vtk

# Lecture du fichier STRUCTURED_GRID (.vtk ASCII)
reader = vtk.vtkStructuredGridReader()
reader.SetFileName("flxmtt_1_1_1.vtk")
reader.Update()

# Écriture en format .vts (VTK XML Structured Grid)
writer = vtk.vtkXMLStructuredGridWriter()
writer.SetFileName("flxmtt_1_1_1.vts")
writer.SetInputData(reader.GetOutput())
writer.SetDataModeToBinary()  # Pour ASCII : writer.SetDataModeToAscii()
writer.SetCompressorTypeToZLib()  # Compression zlib
writer.Write()