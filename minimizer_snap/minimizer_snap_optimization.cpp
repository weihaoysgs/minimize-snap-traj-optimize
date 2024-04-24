#include "minimizer_snap_optimization.h"

// define factorial function, input i, output i!
int MinimizerSnapTrajectoryOptimization::Factorial(int x) {
  int fac = 1;
  for (int i = x; i > 0; i--)
    fac = fac * i;
  return fac;
}

Eigen::MatrixXd MinimizerSnapTrajectoryOptimization::PolyQPGeneration(
        const int d_order,           // the order of derivative
        const Eigen::MatrixXd &Path, // waypoints coordinates (3d)
        const Eigen::MatrixXd &Vel,  // boundary velocity
        const Eigen::MatrixXd &Acc,  // boundary acceleration
        const Eigen::VectorXd &Time) // time allocation in each segment
{
  int p_order = 2 * d_order - 1;         // the order of polynomial
  int p_num1d = p_order + 1;             // the number of variables in each segment
  int m = static_cast<int>(Time.size()); // the number of segments

  // position(x,y,z), so we need (3 * p_num1d) coefficients
  Eigen::MatrixXd PolyCoeff = Eigen::MatrixXd::Zero(m, 3 * p_num1d);

  // compute Q matrix
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(p_num1d * m, p_num1d * m);
  for (int k = 0; k < m; k++) {
    Eigen::MatrixXd Q_seg = Eigen::MatrixXd::Zero(p_num1d, p_num1d);
    for (int i = d_order; i < p_num1d; i++) {
      for (int j = d_order; j < p_num1d; j++) {
        int factor = i + j - 2 * d_order + 1;
        // clang-format off
        double i_fac = static_cast<double>(Factorial(i)) / Factorial(i - d_order);
        double j_fac = static_cast<double>(Factorial(j)) / Factorial(j - d_order);
        Q_seg(i, j) = i_fac * j_fac * 1.0 / (factor)*std::pow(Time(k), factor);
        // clang-format on
      }
    }
    Q.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = Q_seg;
  }

  // compute M matrix, mapping optimization variables into derivatives of each
  // order of state quantities. maps polynomial coefficients to derivatives.
  // TODO only support minimize snap, not support minimize jerk
  assert(d_order == 4);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m * p_num1d, m * p_num1d);

  for (int k = 0; k < m; k++) {
    Eigen::MatrixXd M_k = Eigen::MatrixXd::Zero(p_num1d, p_num1d);
    double t = Time(k);
    double t2 = t * t;
    double t3 = t * t * t;
    double t4 = t2 * t2;
    double t5 = t4 * t;
    double t6 = t4 * t2;
    double t7 = t4 * t3;
    // clang-format off
    M_k << 1 ,   0 ,   0 ,   0 ,   0 ,   0 ,   0 ,   0 ,
           0 ,   1 ,   0 ,   0 ,   0 ,   0 ,   0 ,   0 ,
           0 ,   0 ,   2 ,   0 ,   0 ,   0 ,   0 ,   0 ,
           0 ,   0 ,   0 ,   6 ,   0 ,   0 ,   0 ,   0 ,
           1 ,   t ,   t2,   t3,   t4,   t5,   t6,   t7,
           0 ,   1 ,  2*t,  3*t2, 4*t3, 5*t4, 6*t5,  7*t6,
           0 ,   0 ,   2 ,  6*t ,12*t2,20*t3,30*t4, 42*t5,
           0 ,   0 ,   0 ,  6   , 24*t,60*t2,120*t3,210*t4;
    // clang-format on
    M.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = M_k;
  }

  int C_t_rows = 2 * m * d_order;
  int C_t_cols = 2 * m * d_order - (m - 1) * d_order;
  Eigen::MatrixXd C_T = Eigen::MatrixXd::Zero(C_t_rows, C_t_cols);

  // fix start point state
  for (int i = 0; i < d_order; i++) {
    C_T(i, i) = 1;
  }
  // fix end point state
  int start_col = d_order + m - 1;
  for (int i = C_t_rows - d_order; i < C_t_rows; i++) {
    C_T(i, start_col) = 1;
    start_col++;
  }
  // fixed the position of the intermediate point, implying continuity
  // constraints on the position
  for (int i = 1; i < m; i++) {
    int row = i * 2 * d_order - d_order;
    int col = i + d_order - 1;
    C_T(row, col) = 1;
    row = row + d_order;
    C_T(row, col) = 1;
  }
  // other state constraints at intermediate points, implicit continuity
  // constraints
  for (int i = 1; i < m; i++) {
    int row = i * 2 * d_order - d_order + 1;
    int col = 2 * d_order + m - 1 + (d_order - 1) * (i - 1);
    C_T.block(row, col, d_order - 1, d_order - 1).setIdentity();
    row = row + d_order;
    C_T.block(row, col, d_order - 1, d_order - 1).setIdentity();
  }

  Eigen::MatrixXd C = C_T.transpose();
  Eigen::MatrixXd M_inv = M.inverse();
  Eigen::MatrixXd M_inv_T = M_inv.transpose();
  Eigen::MatrixXd R = C * M_inv_T * Q * M_inv * C_T;

  // Calculate the corresponding coefficients for the three axes respectively
  int num_dim_F = 2 * d_order + m - 1;
  int num_dim_P = static_cast<int>(C_T.cols()) - num_dim_F;

  Eigen::MatrixXd R_pp = R.bottomRightCorner(num_dim_P, num_dim_P);
  Eigen::MatrixXd R_fp = R.topRightCorner(num_dim_F, num_dim_P);

  for (int dim = 0; dim < 3; dim++) {
    Eigen::VectorXd way_pts = Path.col(dim);
    Eigen::VectorXd d_F = Eigen::VectorXd ::Zero(num_dim_F);
    // TODO Fixed starting point state, only the position is fixed here
    d_F(0) = way_pts(0);
    // Fixed position of midpoint
    for (int k = 0; k < m - 1; k++) {
      d_F(d_order + k) = way_pts(k + 1);
    }
    // TODO Fixed ending point state, only the position is fixed here
    d_F(d_order + m - 1) = way_pts(m);
    Eigen::VectorXd d_P = -1.0 * R_pp.inverse() * R_fp.transpose() * d_F;
    Eigen::VectorXd d(num_dim_F + num_dim_P);
    d << d_F, d_P;
    // M*p=C^T*d
    Eigen::VectorXd poly_coef = M_inv * C_T * d;
    Eigen::MatrixXd poly_coef_1d_t = poly_coef.transpose();
    for (int k = 0; k < m; k++) {
      PolyCoeff.block(k, dim * p_num1d, 1, p_num1d) =
              poly_coef_1d_t.block(0, k * p_num1d, 1, p_num1d);
    }
  }
  return PolyCoeff;
}

Eigen::MatrixXd MinimizerSnapTrajectoryOptimization::PolyQPGeneration(
        const int d_order,           // the order of derivative
        const Eigen::MatrixXd &Path, // waypoints coordinates (3d)
        const Eigen::MatrixXd &Vel,  // boundary velocity
        const Eigen::MatrixXd &Acc,  // boundary acceleration
        const Eigen::VectorXd &Time, // time allocation in each segment
        OsqpEigen::Solver &solver) {
  // polynomial order
  int poly_order = 2 * d_order - 1;
  // polynomial parameter number
  int poly_param_num1d = poly_order + 1;
  // trajectory segment
  int m = static_cast<int>(Time.size());
  Eigen::MatrixXd poly_coefficient =
          Eigen::MatrixXd::Zero(m, 3 * poly_param_num1d);

  // set variable param number
  solver.data()->setNumberOfVariables(m * poly_param_num1d);

  // set constraint number. start and end constraint(d_order * 2), middle point
  // constraint(m - 1), continuous constraint((m - 1) * (d_order))
  solver.data()->setNumberOfConstraints(d_order * 2 + (m - 1) * (d_order + 1));

  Eigen::SparseMatrix<double> hessian;
  hessian.resize(m * poly_param_num1d, m * poly_param_num1d);
  hessian.setZero();

  for (int k = 0; k < m; k++) {
    // clang-format off
    Eigen::MatrixXd Q_k = Eigen::MatrixXd::Zero(poly_param_num1d, poly_param_num1d);
    for (int i = d_order; i < poly_param_num1d; i++) {
      for (int j = d_order; j < poly_param_num1d; j++) {
        int factor = i + j - 2 * d_order + 1;
        double i_fac = static_cast<double>(Factorial(i)) / Factorial(i - d_order);
        double j_fac = static_cast<double>(Factorial(j)) / Factorial(j - d_order);
        double value = i_fac * j_fac * 1.0 / (factor)*std::pow(Time(k), factor);
        Q_k(i, j) = value;
        hessian.insert(k * poly_param_num1d + i, k * poly_param_num1d + j) = value;
      }
    }
    // clang-format on
  }

  if (!solver.data()->setHessianMatrix(hessian)) {
    LOG(FATAL) << "Set Hessian Matrix Failed";
  }

  Eigen::SparseMatrix<double> linear_matrix;
  linear_matrix.resize(2 * d_order + (m - 1) * (d_order + 1),
                       m * poly_param_num1d);

  auto generate_coefficient = [&](int row, int col, double t,
                                  bool one_line = false,
                                  bool reverse = false) -> void {
    int flag = d_order;
    if (one_line) {
      flag = 1;
    }
    Eigen::MatrixXd coefficient(d_order, poly_param_num1d);
    if (d_order == 4) {
      coefficient << 1.0, 1.0 * t, 1.0 * pow(t, 2), 1.0 * pow(t, 3),
              1.0 * pow(t, 4), 1.0 * pow(t, 5), 1.0 * pow(t, 6), 1.0 * pow(t, 7),
              0.0, 1.0, 2.0 * t, 3.0 * pow(t, 2), 4.0 * pow(t, 3), 5.0 * pow(t, 4),
              6.0 * pow(t, 5), 7.0 * pow(t, 6), 0.0, 0.0, 2.0, 6.0 * t,
              12.0 * pow(t, 2), 20.0 * pow(t, 3), 30.0 * pow(t, 4),
              42.0 * pow(t, 5), 0.0, 0.0, 0.0, 6.0, 24.0 * t, 60.0 * pow(t, 2),
              120.0 * pow(t, 3), 210.0 * pow(t, 4);
    } else if (d_order == 3) {
      coefficient << 1.0, 1.0 * t, 1.0 * pow(t, 2), 1.0 * pow(t, 3),
              1.0 * pow(t, 4), 1.0 * pow(t, 5), 0.0, 1.0, 2.0 * t, 3.0 * pow(t, 2),
              4.0 * pow(t, 3), 5.0 * pow(t, 4), 0.0, 0.0, 2.0, 6.0 * t,
              12.0 * pow(t, 2), 20.0 * pow(t, 3);
    } else {
      LOG(FATAL)
              << "Only Support Minimizer Jerk/Snap, Please Select Right d_order";
    }
    if (reverse) {
      coefficient = coefficient * (-1.0);
    }
    for (int i = 0; i < d_order && i < flag; ++i) {
      for (int j = 0; j < poly_param_num1d; ++j) {
        linear_matrix.insert(row + i, col + j) = coefficient(i, j);
      }
    }
  };

  int row = 0, col = 0;
  // start point constraint
  generate_coefficient(row, col, 0);
  // end point constraint
  row += d_order;
  col = (m - 1) * poly_param_num1d;
  generate_coefficient(row, col, Time(m - 1));
  row += d_order;
  // middle point position constraint
  for (int k = 0; k < m - 1; k++) {
    generate_coefficient(row + k, k * poly_param_num1d, Time(k), true);
  }
  row += (m - 1);
  // continuous constraint
  for (int k = 0; k < m - 1; k++) {
    generate_coefficient(row, k * poly_param_num1d, Time(k), false, false);
    generate_coefficient(row, (k + 1) * poly_param_num1d, 0, false, true);
    row += d_order;
  }
  if (!solver.data()->setLinearConstraintsMatrix(linear_matrix)) {
    std::cerr << "setLinearConstraintsMatrix Failed\n";
    exit(-1);
  }

  Eigen::VectorXd gradient(poly_param_num1d * m);
  gradient.setZero();

  if (!solver.data()->setGradient(gradient)) {
    LOG(FATAL) << "setGradient Failed";
  }

  Eigen::VectorXd lowbound =
          Eigen::VectorXd::Zero(d_order * 2 + (m - 1) * (d_order + 1));
  Eigen::VectorXd upbound =
          Eigen::VectorXd::Zero(d_order * 2 + (m - 1) * (d_order + 1));

  solver.data()->setLowerBound(lowbound);
  solver.data()->setUpperBound(upbound);

  if (!solver.isInitialized()) {
    solver.initSolver();
  }

  for (int dim = 0; dim < 3; dim++) {
    Eigen::VectorXd way_pts = Path.col(dim);
    lowbound(0) = way_pts(0);
    upbound(0) = way_pts(0);

    lowbound(d_order) = way_pts(m);
    upbound(d_order) = way_pts(m);

    // fix middle point position
    for (int i = 0; i < m - 1; i++) {
      lowbound(2 * d_order + i) = way_pts(i + 1);
      upbound(2 * d_order + i) = way_pts(i + 1);
    }

    // update bounds
    solver.updateBounds(lowbound, upbound);
    // solve QP problem
    solver.solveProblem();

    Eigen::VectorXd poly_coef_1d = solver.getSolution();
    Eigen::MatrixXd poly_coef_1d_t = poly_coef_1d.transpose();

    for (int k = 0; k < m; k++) {
      poly_coefficient.block(k, dim * poly_param_num1d, 1, poly_param_num1d) =
              poly_coef_1d_t.block(0, k * poly_param_num1d, 1, poly_param_num1d);
    }
  }

  solver.data()->clearHessianMatrix();
  solver.data()->clearLinearConstraintsMatrix();
  solver.clearSolverVariables();
  solver.clearSolver();

  return poly_coefficient;
}

Eigen::VectorXd MinimizerSnapTrajectoryOptimization::TimeAllocation(Eigen::MatrixXd Path, double Vel, double Acc) {
  Eigen::VectorXd time(Path.rows() - 1);
  double reach_vel_time_cost = Vel / Acc;
  double distance_acc =
          0.5 * Acc * reach_vel_time_cost * reach_vel_time_cost * 2.0;

  for (Eigen::Index i = 0; i < Path.rows() - 1; i++) {
    double dis = (Path.row(i) - Path.row(i + 1)).norm();
    if (dis <= distance_acc) {
      time(i) = std::sqrt(dis / Acc);
    } else {
      time(i) = 2.0 * reach_vel_time_cost + (dis - distance_acc) / Vel;
    }
  }

  return time;
}

Eigen::Vector3d MinimizerSnapTrajectoryOptimization::getPosPoly(Eigen::MatrixXd polyCoeff, int k, double t, int poly_num1D) {
  Eigen::Vector3d ret;

  for (int dim = 0; dim < 3; dim++) {
    Eigen::VectorXd coeff = (polyCoeff.row(k)).segment(dim * poly_num1D, poly_num1D);
    Eigen::VectorXd time = Eigen::VectorXd::Zero(poly_num1D);

    for (int j = 0; j < poly_num1D; j++)
      if (j == 0)
        time(j) = 1.0;
      else
        time(j) = pow(t, j);

    ret(dim) = coeff.dot(time);
  }

  return ret;
}
