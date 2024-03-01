#ifndef EIGEN_TYPES_H
#define EIGEN_TYPES_H

#include "../LIBS/Eigen/Core"

//typedef Eigen::VectorXd vec;
//typedef Eigen::MatrixXd mat;
//typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> phi_type;
//typedef Eigen::Matrix<double,Eigen::Dynamic,1> t_type;
//typedef Eigen::Matrix3d stress_type;
//typedef Eigen::Matrix<double, Eigen::Dynamic, 6> sd_type;   // sigma_ij order: 11,12,13,22,23,33

/** Atomic Energy Descriptor (AED) type
 *
 * Dimensions:
 * [descriptor size]
 *
 * dim+1 when bias is used
 *
 */

/** Force Descriptor (FD) type
 *
 * Dimensions:
 * [descriptor size, 3]
 *
 */
//typedef Eigen::Matrix<double, Eigen::Dynamic, 3> fd_type;

/** Stress Descriptor (SD) type
 *
 * Dimensions:
 * [descriptor size,6]
 *
 */

/** Density array to be used by many-body descriptors...
 *
 * Dimensions:
 * [atom i][n]
 *
 * where n is specific to given descriptor
 *
 */
//typedef Eigen::Ref<Eigen::VectorXd> rho_type&;
//typedef Eigen::MatrixXd rhos_type;

#endif
