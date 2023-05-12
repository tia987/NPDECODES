/**
 * @file upwindexercise_test.cc
 * @brief NPDE homework UpwindExercise code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "../upwindexercise.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace UpwindExercise::test {

TEST(UpwindExercise, dummyFunction) {
  double x = 0.0;
  int n = 0;

  Eigen::Vector2d v = UpwindExercise::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace UpwindExercise::test
