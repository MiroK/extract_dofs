#include "test.h"     
#include "use_meshfunction.h"
#include "use_indexsets.h"
#include <iostream>

int main()
{
  std::vector<std::string> motion_types;
  motion_types.push_back("rotation");
  motion_types.push_back("translation");

  size_t num_cells[3] = {75, 125, 200};
  size_t widths[4] = {0, 1, 2, 3};

  for(int i = 0; i < motion_types.size(); i++)
  {
    std::cout << motion_types[i] << std::endl;
    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        timing_test(&get_dofs1, &get_dofs0, num_cells[j], widths[k],\
                    motion_types[i]);
      }
    }
  }

  return 0;
}
