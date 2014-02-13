"""
Given a scalar function in space V(mesh) you want to identify dofs that lie
in cells that are part of the band of width 2*w around the contour_value.

1) Find the cells that have the contour_value
2) Look through mesh marking the neighbors of contour_cell. Repeat for depth.
3) Extract dofs all marked cells
"""

from dolfin import *
from numpy import zeros, uintp

def get_dofs0(f, contour_tol, width):
  master_timer = Timer("Getting dofs0")
  master_timer.start()

  mesh = f.function_space().mesh()
  mesh_f = CellFunction("size_t", mesh, 0)

  # get cell with countour
  timer = Timer("First search0")
  timer.start()
  for cell in cells(mesh):
    M = cell.midpoint()
    if contour_tol[0] < f(M) < contour_tol[1]:
      mesh_f[cell] = 1
  timer.stop()
  
  # create a band around cells with countour
  D = mesh.topology().dim()
  N = mesh.num_cells()

  timer = Timer("Connectivity0")
  timer.start()
  mesh.init(D, D)
  timer.stop()

  timer = Timer("Band build0")
  timer.start()
  for depth in range(width):
    new_cells = zeros(N, dtype=uintp) 
    for cell in cells(mesh):
      if mesh_f[cell]:
        neighbors = cell.entities(D)
        for neighbor in neighbors:
          new_cells[int(neighbor)] = 1
   
    mesh_f.set_values(new_cells)
  timer.stop()

  # extract dofs
  timer = Timer("Extract dofs0")
  timer.start()
  dofmap = f.function_space().dofmap()
  dofs = []
  for cell in cells(mesh):
    if mesh_f[cell]:
      cell_dofs = dofmap.cell_dofs(cell.index())
      dofs.extend(cell_dofs)
  timer.stop()

  timer = Timer("Make unique0")
  timer.start()
  dofs = set(dofs) # remove the duplicit dofs
  timer.stop()

  master_timer.stop()

  return dofs
