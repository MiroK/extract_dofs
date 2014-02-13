#ifndef _TEST_H_
#define _TEST_H_

#include <dolfin/common/Array.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <string>
#include <vector>

using namespace dolfin;

class MovingCircle : public Expression
{
public:
  // constructor; circle has radius _r ans is at [_cx, _cy] and will be
  // shifted by _dt everyime update is called
  MovingCircle(double _cx, double _cy, double _r, double _dt) : Expression(),
                                         cx(_cx), cy(_cy), r(_r), dt(_dt) { }
  
  // move the circle; as traslation with velocity (1, 0) or as rotation
  // with velocity (-y, x)
  void update(std::string motion_type);

  // evalaute at point x
  void eval(Array<double>& values, const Array<double>& x) const;

private:
  double cx, cy, r, dt;
};

//-----------------------------------------------------------------------------

// print out timing of get_dofs* methods
void timing_test(std::set<la_index> (*method0)(const Function&,\
                                              const std::vector<double>&,\
                                              const unsigned int),\
                 std::set<la_index> (*method1)(const Function&,\
                                              const std::vector<double>&,\
                                              const unsigned int),\
                 size_t num_cells, size_t width, std::string motion_type);

#endif // _TEST_H_
