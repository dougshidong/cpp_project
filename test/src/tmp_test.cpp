#include "project/tmp.hpp"

#include <gtest/gtest.h>

TEST(TmpAddTest, CheckValues)
{
  ASSERT_EQ(tmp::add(1, 2), 3);
}

TEST(ConvectivePhysicsTest, TmpAddTest_CheckValues_Test) {
  // Tests ConvectivePhysics
  EulerParams const params =
      EulerParams::Builder().set_gamma(1.4).set_ref_temperature(1.0).set_ref_density(1.0).set_ref_pressure(1.0).build();
  constexpr int dim = 1;
  constexpr int nstate = dim+2;
  Euler<dim, nstate, double> euler(params);
  euler.print();
  std::array<double, nstate> const state = { 1.0, 2.0, 3.0 };
  std::array<std::array<double, nstate>, dim> flux = euler.compute_flux(state);
  for (size_t i = 0; i < dim; i++) {
    std::cout << flux[i][0] << std::endl;
  }
  // Tests ConvectivePhysics
  EXPECT_DOUBLE_EQ(flux[0][0], 2.0);
  EXPECT_DOUBLE_EQ(flux[0][1], 4.4);
  EXPECT_DOUBLE_EQ(flux[0][2], 19.8);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
