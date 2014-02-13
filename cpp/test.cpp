#include "test.h"
#include "Poisson.h"
#include <dolfin/generation/RectangleMesh.h>
#include <dolfin/log/log.h>
#include "dolfin/common/timing.h"

void MovingCircle::update(std::string motion_type)
{
  if(motion_type == "translation")
  {
    cx += dt;
  }
  else if(motion_type == "rotation")
  {
    double _cx = cx;
    cx -= cy*dt;
    cy += _cx*dt;
  }
}

//-----------------------------------------------------------------------------

void MovingCircle::eval(Array<double>& values, const Array<double>& x) const
{
  values[0] = r - sqrt(pow(x[0] - cx, 2) + pow(x[1] - cy, 2));
}

//-----------------------------------------------------------------------------

void timing_test(std::set<la_index> (*method0)(const Function&,\
                                              const std::vector<double>&,\
                                              const unsigned int),\
                 std::set<la_index> (*method1)(const Function&,\
                                              const std::vector<double>&,\
                                              const unsigned int),\
                 size_t num_cells, size_t width, std::string motion_type)
{
  RectangleMesh mesh(-2, -2, 2, 2, num_cells, num_cells);
  mesh.init(2, 2);
  Poisson::FunctionSpace V(mesh);

  MovingCircle circle(-1, 0, 0.5, 0.1);
  double _contour_tol[2] = {-0.1, 0.1};
  std::vector<double> contour_tol(_contour_tol, _contour_tol + 2);  

  Function u(V);

  std::set<la_index> dofs0;
  std::set<la_index> dofs1;
  for(int i = 0; i < 10; i++)
  {
    circle.update(motion_type);
    u.interpolate(circle);

    dofs0 = (*method0)(u, contour_tol, width);
    dofs1 = (*method1)(u, contour_tol, width);
  }

  // check the difference of sets 
  std::set<la_index> difference;
  std::set_difference(dofs0.begin(), dofs0.end(), dofs1.begin(), dofs1.end(),
                      std::inserter(difference, difference.end()));
  info("Dofs differ by set of measure %d", difference.size());
  
  // check the performance
  info("width %d, n_cells %d cells, time %g vs %f in dofs of which",\
  2*width, mesh.num_cells(), timing("Getting dofs0"), timing("Getting dofs1"));
  
  info("\t\t first search %g vs %g",\
  timing("First search0"), timing("First search1"));

  info("\t\t band build %g vs %g",\
  timing("Band build0"), timing("Band build1"));

  info("\t\t connectivity %g vs %g",\
  timing("Connectivity0"), timing("Connectivity1"));

  info("\t\t extract dofs %g vs %g",\
  timing("Extract dofs0"), timing("Extract dofs1"));
}

