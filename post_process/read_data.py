import numpy as np
import vtk
from vtk.util import numpy_support as VN


class read_data:

  def __init__(self):
      pass

  def read_vtk(self, fname):

      reader = vtk.vtkDataSetReader()

      reader.SetFileName(fname)
      reader.ReadAllScalarsOn()
      reader.Update()

      ### Read solution fields
      data = reader.GetOutput()
      dim = data.GetDimensions()
      vec = list(reversed(list(dim)))[1:] # Ny x Nx
      Ny = vec[0]
      Nx = vec[1]

      # Solution fields
      ux = VN.vtk_to_numpy(data.GetPointData().GetArray('ux'))
      uy = VN.vtk_to_numpy(data.GetPointData().GetArray('uy'))
      T = VN.vtk_to_numpy(data.GetPointData().GetArray('Temperature'))

      # Reshape them
      ux = ux.reshape(vec)
      uy = uy.reshape(vec)
      T = T.reshape(vec)

      # Get mesh points
      x = np.zeros(data.GetNumberOfPoints())
      y = np.zeros(data.GetNumberOfPoints())
      z = np.zeros(data.GetNumberOfPoints())
      
      for i in range(data.GetNumberOfPoints()):
              x[i],y[i],z[i] = data.GetPoint(i)
      
      X = x.reshape(vec)
      Y = y.reshape(vec)
      Z = z.reshape(vec)

      return X, Y, ux, uy, T
