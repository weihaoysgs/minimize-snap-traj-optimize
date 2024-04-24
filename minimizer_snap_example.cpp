#include "minimizer_snap/matplotlibcpp.hpp"
#include "minimizer_snap/minimizer_snap_optimization.h"
#include <glog/logging.h>

OsqpEigen::Solver osqp_solver;

namespace plt = matplotlibcpp;

std::vector<Eigen::Vector3d> TrajGeneration(const Eigen::MatrixXd &path) {
  int d_order = 4;
  Eigen::MatrixXd vel = Eigen::MatrixXd::Zero(2, 3);
  Eigen::MatrixXd acc = Eigen::MatrixXd::Zero(2, 3);
  Eigen::VectorXd poly_time = MinimizerSnapTrajectoryOptimization::TimeAllocation(path, 1.0, 1.0);
  // you can add the last parameter for osqp solver to using OSQP iterative solve the QP problem
  Eigen::MatrixXd poly_coefficient = MinimizerSnapTrajectoryOptimization::PolyQPGeneration(
          d_order, path, vel, acc, poly_time/*,osqp_solver*/);
  std::cout << "\n" << poly_coefficient;
  std::vector<Eigen::Vector3d> smooth_traj_path;
  for (int i = 0; i < poly_time.size(); i++) {
    for (double t = 0.0; t < poly_time(i);) {
      Eigen::Vector3d pos = MinimizerSnapTrajectoryOptimization::getPosPoly(poly_coefficient, i, t, 2 * d_order);
      t += 0.01;
      smooth_traj_path.emplace_back(pos);
    }
  }
  return smooth_traj_path;
}

int main() {
  google::InitGoogleLogging("minimizer_snap_example");
  FLAGS_colorlogtostderr = true;
  FLAGS_stderrthreshold = google::INFO;
  Eigen::Matrix<double, 7, 3> waypoints;
  // clang-format off
  waypoints <<         0,         0,      0,
                 1.20292,   3.05196,   1.32,
               -0.692698,   2.74029,   0.88,
               -0.979247, -0.258448,   1.46,
                0.366166,  -1.73207,      2,
                 3.47916,  -1.59115,   3.62,
                 3.78978,  0.978379,   0.56;
  // clang-format on

  // solver the minimizer snap problem
  std::vector<Eigen::Vector3d> smooth_traj_path = TrajGeneration(waypoints);

  std::vector<double> x, smooth_x, y, smooth_y, z, smooth_z;
  std::map<std::string, std::string> keywords, smooth_keywords;
  keywords.insert(std::pair<std::string, std::string>("label", "Way Points"));
  smooth_keywords.insert(std::pair<std::string, std::string>("label", "Minimizer Snap Traj"));
  smooth_keywords.insert(std::pair<std::string, std::string>("color", "red"));
  for (auto i = 0; i < waypoints.rows(); i++) {
    x.push_back(waypoints.row(i)(0));
    y.push_back(waypoints.row(i)(1));
    z.push_back(waypoints.row(i)(2));
  }
  for (auto i = 0; i < smooth_traj_path.size(); i++) {
    smooth_x.push_back(smooth_traj_path[i].x());
    smooth_y.push_back(smooth_traj_path[i].y());
    smooth_z.push_back(smooth_traj_path[i].z());
  }

  plt::plot3(x, y, z, smooth_x, smooth_y, smooth_z, keywords, smooth_keywords);
  plt::xlabel("x label");
  plt::ylabel("y label");
  plt::set_zlabel("z label");
  plt::legend();
  plt::show();
  plt::detail::_interpreter::kill();
  return 0;
}


