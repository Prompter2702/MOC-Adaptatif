import vtk

# Lecture du fichier STRUCTURED_GRID (.vtk ASCII)
reader = vtk.vtkStructuredGridReader()
reader.SetFileName("flx_vol.vtk")
reader.Update()

# Écriture en format .vts (VTK XML Structured Grid)
writer = vtk.vtkXMLStructuredGridWriter()
writer.SetFileName("flx_vol.vts")
writer.SetInputData(reader.GetOutput())
writer.SetDataModeToBinary()  # Pour ASCII : writer.SetDataModeToAscii()
writer.SetCompressorTypeToZLib()  # Compression zlib
writer.Write()