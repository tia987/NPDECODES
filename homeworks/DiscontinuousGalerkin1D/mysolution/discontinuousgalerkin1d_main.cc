/**
 * @file discontinuousgalerkin1d_main.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "discontinuousgalerkin1d.h"

int main() {
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file;
  file.open("solution.csv");
  for(double T = 0.; T <= 5; T+=1){
      DiscontinuousGalerkin1D::Solution solution = DiscontinuousGalerkin1D::solveTrafficFlow(T);


      //====================
      // Your code goes here
      // Use std::ofstream to write the solution to
      // the file "solution.csv". To plot this file
      // you may uncomment the following line:
      
      if(T < 0.01) file << solution.x_.transpose().format(CSVFormat) << std::endl;
      file << solution.u_.transpose().format(CSVFormat) << std::endl;
  }
  file.close();
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_solution.py "
  CURRENT_BINARY_DIR "/solution.csv " CURRENT_SOURCE_DIR "/solution.jpg");
  //====================

  return 0;
}
