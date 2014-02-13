#ifndef _USE_MESHFUNCTION_H_
#define _USE_MESHFUNCTION_H_

#include <set>
#include <vector>
#include <dolfin/common/types.h>
#include <dolfin/function/Function.h>

using namespace dolfin;

// extract dofs of the scalar function that lie in cells in band of width around
// cells where F takes values between contour_tol[0] and contour_tol[1]
std::set<la_index> get_dofs0(const Function& F,\
                             const std::vector<double>& contour_tol,\
                             const unsigned int width);

#endif // _USE_MESHFUNCTION_H_
