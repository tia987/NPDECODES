/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include "transformationofgalerkinmatrices.h"

using namespace TransformationOfGalerkinMatrices;

int main() {
  std::cout << "NPDE Homework problem: TransformationOfGalerkinMatrices"
            << std::endl;

  std::vector<Eigen::Triplet<double>> A;
  A.push_back(Eigen::Triplet<double>(0, 0, 1));
  A.push_back(Eigen::Triplet<double>(1, 1, 1));

  std::vector<Eigen::Triplet<double>> A_tilde = transformCOOmatrix(A);

  std::cout << "Run the tests to check your code!" << std::endl;

  return 0;
}
