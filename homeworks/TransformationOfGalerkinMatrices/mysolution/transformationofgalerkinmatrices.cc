/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "transformationofgalerkinmatrices.h"

#include <cassert>
#include <iostream>

namespace TransformationOfGalerkinMatrices {

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Triplet<double>> transformCOOmatrix(
    const std::vector<Eigen::Triplet<double>> &A) {
  std::vector<Eigen::Triplet<double>> A_t{};  // return value

  // First step: find the size of the matrix by searching the maximal
  // indices. Depends on the assumption that no zero rows/columns occur.
  int rows_max_idx = 0, cols_max_idx = 0;
  for (const Eigen::Triplet<double> &triplet : A) {
    rows_max_idx =
        (triplet.row() > rows_max_idx) ? triplet.row() : rows_max_idx;
    cols_max_idx =
        (triplet.col() > cols_max_idx) ? triplet.col() : cols_max_idx;
  }
  int n_rows = rows_max_idx + 1;
  int n_cols = cols_max_idx + 1;

  // Make sure we deal with a square matrix
  assert(n_rows == n_cols);
  // The matrix size must have even parity
  assert(n_cols % 2 == 0);

  int N = n_cols;      // Size of (square) matrix
  int M = n_cols / 2;  // Half the size
  //====================
  // Your code goes here
  //====================
  Eigen::SparseMatrix<double> m(N,N);
  m.setFromTriplets(A.begin(), A.end());
  for (size_t i = 1; i < N; i++){
    for (size_t j = 1; j < N; j++){    
      if(i < M && j < M){
        double t = m.coeff(2*i,2*j)+m.coeff(2*i,2*j-1)+m.coeff(2*i-1,2*j)+m.coeff(2*i-1,2*j-1);
        Eigen::Triplet<double> trip(i,j,t);
        A_t.push_back(trip);
      }
      else if(M < i && i < N && M < j && j < N){        
        auto t = m.coeff(2*(i-M)-1,2*(j-M)-1)+m.coeff(2*(i-M)-1,2*(j-M))-m.coeff(2*(i-M),2*(j-M)-1)+m.coeff(2*(i-M),2*(j-M));
        Eigen::Triplet<double> trip(i,j,t);
        A_t.push_back(trip);
      }
      else if(M < i && i < N && 0 < j && j < M){
        auto t = m.coeff(2*(i-M)-1,2*j-1)+m.coeff(2*(i-M)-1,2*j)-m.coeff(2*(i-M),2*j-1)-m.coeff(2*(i-M),2*j);
        Eigen::Triplet<double> trip(i,j,t);
        A_t.push_back(trip);
      }
      else if(0 < i && i < M && M < j && j < N){
        auto t = m.coeff(2*i-1,2*(j-M)-1)-m.coeff(2*i-1,2*(j-M))+m.coeff(2*i,2*(j-M)-1)-m.coeff(2*i,2*(j-M));
        Eigen::Triplet<double> trip(i,j,t);
        A_t.push_back(trip);
      }
    }
  }

  return A_t;
}
/* SAM_LISTING_END_1 */

}  // namespace TransformationOfGalerkinMatrices
