#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  n_aug_ = 7;

  n_x_ = 5;

  is_initialized_ = false;

  lambda_ = 3 - n_aug_;

  time_us_ = 0;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2*n_aug_+1);

  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for(int i = 1; i < (2*n_aug_+1) ; i++)
  {
      weights_(i) = 0.5 /  (lambda_ + n_aug_);
  }

  x_.fill(0.0);

  Xsig_pred_.fill(0.0);

  P_.setIdentity();
  P_(3,3) = 0.5;
  P_(4,4) = 0.5;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

    if(!is_initialized_)
    {
        if(meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            x_ << meas_package.raw_measurements_[0],
                  meas_package.raw_measurements_[1],
                  0,
                  0,
                  0;

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]),
                  meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]),
                  0,
                  0,
                  0;
        }
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }

    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(delta_t);

    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        UpdateRadar(meas_package);
    }
    else
    {
        /* Next Sensor ;) */
    }


}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

    //////////////////////////////////////////////////////////////////
    //GENERATE SIGMA PTS
    //////////////////////////////////////////////////////////////////

    MatrixXd Xsig_gen = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_gen.fill(0.0);

    MatrixXd Q = MatrixXd(2, 2);
    Q <<    std_a_*std_a_,     0,
            0,  std_yawdd_*std_yawdd_;

    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.fill(0.0);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);


    // create augmented mean state
    x_aug.head(n_x_) = x_;

    // create augmented covariance matrix
    P_aug.topLeftCorner(P_.rows(),P_.cols()) = P_;
    P_aug.block(n_x_,n_x_,2,2) = Q;

    // create square root matrix
    MatrixXd A = P_aug.llt().matrixL();

    // create augmented sigma points
    double sqrt_val = sqrt(lambda_+n_aug_);

    Xsig_gen.col(0) = x_aug;

    for(int i =0;i<Xsig_gen.rows();i++)
    {
        Xsig_gen.col(i+1) = x_aug + sqrt_val * A.col(i);
        Xsig_gen.col(i+1+n_aug_) = x_aug - sqrt_val * A.col(i);
    }

    //////////////////////////////////////////////////////////////////
    //PREDICT SIGMA PTS
    //////////////////////////////////////////////////////////////////

    double px, py, v, thai, thai_dot, nu_a, nu_thai;
    VectorXd x_k = VectorXd(5);
    VectorXd f_k = VectorXd(5);
    VectorXd meu_k = VectorXd(5);

    for(int i = 0; i < Xsig_gen.cols() ; i++)
    {
        px = Xsig_gen(0,i);
        py = Xsig_gen(1,i);
        v = Xsig_gen(2,i);
        thai= Xsig_gen(3,i);
        thai_dot = Xsig_gen(4,i);
        nu_a = Xsig_gen(5,i);
        nu_thai = Xsig_gen(6,i);

        if(fabs(thai_dot) > std::numeric_limits<double>::epsilon())
        {
            f_k(0) = (v / thai_dot) * ( sin(thai + thai_dot*delta_t) - sin(thai) );
            f_k(1) = (v / thai_dot) * ( cos(thai) - cos(thai + thai_dot*delta_t) );
            f_k(2) = 0;
            f_k(3) = thai_dot * delta_t;
            f_k(4) = 0;
        }
        else
        {
            f_k(0) = v * cos(thai) * delta_t;
            f_k(1) = v * sin(thai) * delta_t;
            f_k(2) = 0;
            f_k(3) = 0;
            f_k(4) = 0;
        }

        //      std::cout << "f_k = " << std::endl << f_k << std::endl;

        meu_k(0) = 0.5 * delta_t * delta_t * cos(thai) * nu_a;
        meu_k(1) = 0.5 * delta_t * delta_t * sin(thai) * nu_a;
        meu_k(2) = delta_t * nu_a;
        meu_k(3) = 0.5 * delta_t * delta_t * nu_thai;
        meu_k(4) = delta_t * nu_thai;

        x_k = Xsig_gen.block(0,i,n_x_,1);

        Xsig_pred_.col(i) = x_k + f_k + meu_k;
    }

    //////////////////////////////////////////////////////////////////
    //CALCULATE MEAN AND COVARIANCE
    //////////////////////////////////////////////////////////////////

    x_.fill(0.0);
    P_.fill(0.0);

    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        x_ += Xsig_pred_.col(i) * weights_(i);

    }

    // predict state covariance matrix

    VectorXd tmp_vec = VectorXd(n_x_);
    tmp_vec.fill(0.0);
    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        tmp_vec = Xsig_pred_.col(i) - x_;
        while(tmp_vec(3) > M_PI)
        {
            tmp_vec(3) = tmp_vec(3) - 2. * M_PI;
        }
        while(tmp_vec(3) < -M_PI)
        {
            tmp_vec(3) = tmp_vec(3) + 2. * M_PI;
        }
        P_ += weights_(i) * ( tmp_vec * tmp_vec.transpose());
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
     * TODO: Complete this function! Use lidar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the lidar NIS, if desired.
     */

    int n_z = 2;
    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);

    // create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z.fill(0.0);
    z << meas_package.raw_measurements_[0],   // px in m
         meas_package.raw_measurements_[1];   // py in m

    // transform sigma points into measurement space

    double px, py, v, thai;
    for(int i = 0; i < 2*n_aug_+1 ; i++)
    {
        Zsig(0,i) = Xsig_pred_(0,i);
        Zsig(1,i) = Xsig_pred_(1,i);
    }

    // calculate mean predicted measurement

    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        z_pred += weights_(i) * Zsig.col(i);

    }

    // calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z,n_z);
    R.fill(0.0);
    R(0,0) = std_laspx_ * std_laspx_;
    R(1,1) = std_laspy_ * std_laspy_;

    VectorXd tmp_vec = VectorXd(n_z);
    tmp_vec.fill(0.0);
    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        tmp_vec = Zsig.col(i) - z_pred;
        while(tmp_vec(1) > M_PI)
        {
            tmp_vec(1) = tmp_vec(1) - 2. * M_PI;
        }
        while(tmp_vec(1) < -M_PI)
        {
            tmp_vec(1) = tmp_vec(1) + 2. * M_PI;
        }
        S += weights_(i) * ( tmp_vec * tmp_vec.transpose());
    }

    S += R;

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    // create matrix for Kalman Gain
    MatrixXd K_gain = MatrixXd(n_x_, n_z);
    K_gain.fill(0.0);

    // calculate cross correlation matrix
    for(int i=0; i<2*n_aug_+1;i++)
    {
        VectorXd xdiff = VectorXd(n_x_);
        xdiff.fill(0.0);
        VectorXd zdiff = VectorXd(n_z);
        zdiff.fill(0.0);

        xdiff = (Xsig_pred_.col(i) - x_);
        zdiff = (Zsig.col(i) - z_pred);

        while(xdiff(3) > M_PI)
        {
            xdiff(3) = xdiff(3) - 2. * M_PI;
        }
        while(xdiff(3) < -M_PI)
        {
            xdiff(3) = xdiff(3) + 2. * M_PI;
        }


        while(zdiff(1) > M_PI)
        {
            zdiff(1) = zdiff(1) - 2. * M_PI;
        }
        while(zdiff(1) < -M_PI)
        {
            zdiff(1) = zdiff(1) + 2. * M_PI;
        }


        Tc += weights_(i) * xdiff * zdiff.transpose();
    }

    // calculate Kalman gain K;

    K_gain = Tc * S.inverse();

    // update state mean and covariance matrix
    VectorXd zdiff = VectorXd(n_z);
    zdiff.fill(0.0);
    zdiff = z - z_pred;
    while(zdiff(1) > M_PI)
    {
        zdiff(1) = zdiff(1) - 2. * M_PI;
    }
    while(zdiff(1) < -M_PI)
    {
        zdiff(1) = zdiff(1) + 2. * M_PI;
    }

    x_ = x_ + K_gain * zdiff;

    P_ = P_ - K_gain * S * K_gain.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
     * TODO: Complete this function! Use radar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the radar NIS, if desired.
     */

    int n_z = 3;
    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);

    // create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z.fill(0.0);
    z << meas_package.raw_measurements_[0],   // rho in m
            meas_package.raw_measurements_[1],   // phi in rad
            meas_package.raw_measurements_[2];   // rho_dot in m/s

    // transform sigma points into measurement space

    double px, py, v, thai;
    for(int i = 0; i < 2*n_aug_+1 ; i++)
    {
        px = Xsig_pred_(0,i);
        py = Xsig_pred_(1,i);
        v = Xsig_pred_(2,i);
        thai= Xsig_pred_(3,i);

        Zsig(0,i) = sqrt(px*px + py*py);
        Zsig(1,i) = atan2(py,px);
        Zsig(2,i) = (px * cos(thai)*v + py*sin(thai)*v) /  sqrt(px*px + py*py);

    }

    // calculate mean predicted measurement

    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        z_pred += weights_(i) * Zsig.col(i);

    }

    // calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z,n_z);
    R.fill(0.0);
    R(0,0) = std_radr_ * std_radr_;
    R(1,1) = std_radphi_ * std_radphi_;
    R(2,2) = std_radrd_ * std_radrd_;

    VectorXd tmp_vec = VectorXd(n_z);
    tmp_vec.fill(0.0);
    for(int i = 0; i < (2*n_aug_+1) ; i++)
    {
        tmp_vec = Zsig.col(i) - z_pred;
        while(tmp_vec(1) > M_PI)
        {
            tmp_vec(1) = tmp_vec(1) - 2. * M_PI;
        }
        while(tmp_vec(1) < -M_PI)
        {
            tmp_vec(1) = tmp_vec(1) + 2. * M_PI;
        }
        S += weights_(i) * ( tmp_vec * tmp_vec.transpose());
    }

    S += R;

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    // create matrix for Kalman Gain
    MatrixXd K_gain = MatrixXd(n_x_, n_z);
    K_gain.fill(0.0);


    // calculate cross correlation matrix
    for(int i=0; i<2*n_aug_+1;i++)
    {
        VectorXd xdiff = VectorXd(n_x_);
        VectorXd zdiff = VectorXd(n_z);

        xdiff = (Xsig_pred_.col(i) - x_);
        zdiff = (Zsig.col(i) - z_pred);

        while(xdiff(3) > M_PI)
        {
            xdiff(3) = xdiff(3) - 2. * M_PI;
        }
        while(xdiff(3) < -M_PI)
        {
            xdiff(3) = xdiff(3) + 2. * M_PI;
        }


        while(zdiff(1) > M_PI)
        {
            zdiff(1) = zdiff(1) - 2. * M_PI;
        }
        while(zdiff(1) < -M_PI)
        {
            zdiff(1) = zdiff(1) + 2. * M_PI;
        }

        Tc += weights_(i) * xdiff * zdiff.transpose();
    }

    // calculate Kalman gain K;

    K_gain = Tc * S.inverse();

    // update state mean and covariance matrix

    VectorXd zdiff = VectorXd(n_z);
    zdiff.fill(0.0);
    zdiff = z - z_pred;
    while(zdiff(1) > M_PI)
    {
        zdiff(1) = zdiff(1) - 2. * M_PI;
    }
    while(zdiff(1) < -M_PI)
    {
        zdiff(1) = zdiff(1) + 2. * M_PI;
    }

    x_ = x_ + K_gain * zdiff;

    P_ = P_ - K_gain * S * K_gain.transpose();
}