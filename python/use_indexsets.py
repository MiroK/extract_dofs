"""
Given a scalar function in space V(mesh) you want to identify dofs that lie
in cells that are part of the band of width 2*w around the contour_value.

1) Find the cells that have the contour_value
2) Look through mesh marking the neighbors of contour_cell. Repeat for depth.
3) Extract dofs all marked cells
"""

from dolfin import *

def get_dofs1(f, contour_tol, width):
  master_timer = Timer("Getting dofs1")
  master_timer.start()

  mesh = f.function_space().mesh()
  indices = []

  # get cell with countour
  timer = Timer("First search1")
  timer.start()
  for cell in cells(mesh):
    M = cell.midpoint()
    if contour_tol[0] < f(M) < contour_tol[1]:
      indices.append(cell.index()) # at this point indices contains unique elems
  timer.stop()
  
  # create a band around cells with countour
  timer = Timer("Connectivity1")
  timer.start()
  D = mesh.topology().dim()
  mesh.init(D, D)
  timer.stop()

  timer = Timer("Band build1")
  timer.start()
  start = 0
  stop = len(indices)
  for depth in range(width):
    new_indices = []
    for cell in indices[start:stop]:
      new_indices.extend(Cell(mesh, cell).entities(D))
    
    # make unique and really new
    new_indices = list(set(new_indices) - set(indices))
    indices.extend(new_indices)

    start = stop
    stop += len(new_indices)
  timer.stop()

  # extract dofs
  timer = Timer("Extract dofs1")
  timer.start()
  dofmap = f.function_space().dofmap()
  dofs = []
  for cell in indices:
    cell_dofs = dofmap.cell_dofs(int(cell))
    dofs.extend(cell_dofs)
  timer.stop()

  timer = Timer("Make unique1")
  timer.start()
  dofs = set(dofs)             # remove the duplicit dofs
  timer.stop()

  master_timer.stop()  
  
  return dofs

