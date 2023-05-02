#include <iostream>
#include "Chebyshev.h"
#include <ctime>

int main()
{
  Chebyshev* tmp = new Test_function(100, 100, 10);
  size_t time_start = clock();
  tmp->run();
  size_t time_end = clock();
  std::cout << tmp->inf_discrepancy() << " time work: " << (time_end - time_start)/1000.0;
  return 0;
}