#ifndef TRAJECTORY_GENERATOR_WAYPOINT_H_
#define TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <OsqpEigen/OsqpEigen.h>
#include <glog/logging.h>
#include <iostream>
#include <vector>

class MinimizerSnapTrajectoryOptimization {
public:
  MinimizerSnapTrajectoryOptimization() = default;
  ~MinimizerSnapTrajectoryOptimization() = default;

  static Eigen::VectorXd TimeAllocation(Eigen::MatrixXd Path, double _Vel, double _Acc);
  static Eigen::Vector3d getPosPoly(Eigen::MatrixXd polyCoeff, int k, double t, int poly_num1D);

  static Eigen::MatrixXd PolyQPGeneration(int order,
                                          const Eigen::MatrixXd &Path,
                                          const Eigen::MatrixXd &Vel,
                                          const Eigen::MatrixXd &Acc,
                                          const Eigen::VectorXd &Time);

  static Eigen::MatrixXd
  PolyQPGeneration(int d_order, const Eigen::MatrixXd &Path,
                   const Eigen::MatrixXd &Vel, const Eigen::MatrixXd &Acc,
                   const Eigen::VectorXd &Time, OsqpEigen::Solver &solver);
  static int Factorial(int x);
};

#endif
