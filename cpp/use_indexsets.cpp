#include "use_indexsets.h"
#include <dolfin/mesh/Cell.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/common/Timer.h>
#include <boost/range/algorithm/set_algorithm.hpp>

std::set<la_index> get_dofs1(const Function& F,\
                             const std::vector<double>& contour_tol,\
                             const unsigned int width)
{
  Timer master_timer("Getting dofs1");
  master_timer.start();

  assert(F.value_rank() == 0);                    // only scalar valued funcs
  assert(contour_tol.size() == 2);                // check bounds
  assert(contour_tol[0] < contour_tol[1]);

  boost::shared_ptr<const Mesh> mesh(F.function_space()->mesh()); 
  unsigned int gdim = mesh->geometry().dim();
  unsigned int tdim = mesh->topology().dim();
  std::vector<size_t> indices;

  Timer search_timer("First search1");
  search_timer.start();
  // find cells with F(midpoint) within tolerance
  Array<double> values(0);
  
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    Array<double> x(gdim, cell->midpoint().coordinates());
    F.eval(values, x);
    if((values[0] > contour_tol[0]) and (values[0] < contour_tol[1]))
    {
      indices.push_back(cell->index()); 
    }
  }
  search_timer.stop();

  // loop over already marked cells and mark their neigboring cells
  Timer connectivity_timer("Connectivity1");
  connectivity_timer.start();
  mesh->init(tdim, tdim);
  unsigned int N = mesh->num_cells();
  connectivity_timer.stop();
 
  Timer band_timer("Band build1");
  band_timer.start();

  unsigned int start = 0;
  unsigned int len = indices.size() - start;
  unsigned int stop = start + len;

  std::vector<size_t> new_indices;
  std::set<size_t> _new_indices, _indices, unique_new_indices;
  
  for(unsigned int depth = 0; depth < width; depth++)
  {
    new_indices.clear();
    for(unsigned int  cell = start; cell < stop; cell++)
    {
      Cell _cell(*mesh, indices[cell]);
      const unsigned int* neighbors = _cell.entities(tdim);
      const unsigned int n_neighbors = _cell.num_entities(tdim);
      new_indices.insert(new_indices.end(), neighbors, neighbors + n_neighbors);
    }

    // remove those new indices that are already in indices and add the unique
    // ones
    _new_indices.clear(); _indices.clear(); unique_new_indices.clear();

    _new_indices.insert(new_indices.begin(), new_indices.end());
    _indices.insert(indices.begin(), indices.end());
    boost::set_difference(_new_indices, _indices, std::back_inserter(indices));
    
    start = stop;
    len = indices.size() - start;
    stop = start + len;
  }
  band_timer.stop();

  Timer extract_timer("Extract dofs1");
  extract_timer.start();
  // extract dofs of marked cells
  std::set<la_index> dofs;

  boost::shared_ptr<const GenericDofMap> dofmap(F.function_space()->dofmap());
  for(std::vector<size_t>::const_iterator cell = indices.begin();\
      cell != indices.end(); cell++)
  {
    std::vector<la_index> cell_dofs = dofmap->cell_dofs(*cell);
    dofs.insert(cell_dofs.begin(), cell_dofs.end());
  }
  extract_timer.stop();

  master_timer.stop();
  return dofs;
}
