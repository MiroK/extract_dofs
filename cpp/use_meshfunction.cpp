#include "use_meshfunction.h"
#include <dolfin/mesh/Cell.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/common/Timer.h>

std::set<la_index> get_dofs0(const Function& F,\
                             const std::vector<double>& contour_tol,\
                             const unsigned int width)
{
  Timer master_timer("Getting dofs0");
  master_timer.start();

  assert(F.value_rank() == 0);                    // only scalar valued funcs
  assert(contour_tol.size() == 2);                // check bounds
  assert(contour_tol[0] < contour_tol[1]);

  boost::shared_ptr<const Mesh> mesh(F.function_space()->mesh()); 
  CellFunction<bool> cell_f(*mesh, false);
  unsigned int gdim = mesh->geometry().dim();
  unsigned int tdim = mesh->topology().dim();

  Timer search_timer("First search0");
  search_timer.start();
  // find cells with F(midpoint) within tolerance
  Array<double> values(0);
  
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    Array<double> x(gdim, cell->midpoint().coordinates());
    F.eval(values, x);
    if((values[0] > contour_tol[0]) and (values[0] < contour_tol[1]))
    {
      cell_f[*cell] = true; 
    }
  }
  search_timer.stop();

  // loop over already marked cells and mark their neigboring cells
  Timer connectivity_timer("Connectivity0");
  connectivity_timer.start();
  mesh->init(tdim, tdim);
  unsigned int N = mesh->num_cells();
  connectivity_timer.stop();
 
  Timer band_timer("Band build0");
  band_timer.start();
  for(unsigned int depth = 0; depth < width; depth++)
  {
    std::vector<bool> new_cells(N, false);
    for(CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      if(cell_f[*cell])
      {
        const unsigned int* neighbors = cell->entities(tdim);
        const unsigned int n_neighbors = cell->num_entities(tdim);
        for(unsigned int i = 0; i < n_neighbors; i++)
        { 
          new_cells[neighbors[i]] = true;
        }
      }
    }
    cell_f.set_values(new_cells);
  }
  band_timer.stop();

  Timer extract_timer("Extract dofs0");
  extract_timer.start();
  // extract dofs of marked cells
  std::set<la_index> dofs;

  boost::shared_ptr<const GenericDofMap> dofmap(F.function_space()->dofmap());
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    if(cell_f[*cell])
    {
      std::vector<la_index> cell_dofs = dofmap->cell_dofs(cell->index());
      dofs.insert(cell_dofs.begin(), cell_dofs.end());
    }
  }
  extract_timer.stop();

  master_timer.stop();
  return dofs;
}
