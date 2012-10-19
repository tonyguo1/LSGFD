/*
 * MHDSolver.h
 *
 *  Created on: Oct 19, 2012
 *      Author: tongfei
 */

#ifndef MHDSOLVER_H_
#define MHDSOLVER_H_

#include "EllipticSolverBase.h"
#include <vector>
#include "solver_lapack_cf.h"

namespace std {

class MHD_Solver: public std::Elliptic_Solver_Base {
public:
	MHD_Solver(DATA *data);
	virtual ~MHD_Solver();
	//! The method to solve elliptic problem
	virtual void Solve();
	//! The method for setting Matrix
	virtual void SetMatrix(const vector<vector<int> > &neighbour_list,
				           const vector<vector<double> > &coefficient_laplacian,
			               const vector<vector<double> > &coefficient_dudx,
			               const vector<vector<double> > &coefficient_dudy,
			               const vector<vector<double> > &coefficient_dudz,
			               const vector<vector<double> > &coefficient_laplacian_boundary,
			               const vector<vector<double> > &coefficient_dudx_boundary,
			               const vector<vector<double> > &coefficient_dudy_boundary,
			               const vector<vector<double> > &coefficient_dudz_boundary,
			               const vector<double> &normal_x,
			               const vector<double> &normal_y,
			               const vector<double> &normal_z);
	//! The method to initialize UxB
	void Initializae_UxB();
	//! The method to set Phi
	void Set_Phi(vector<double> &phi, double *x);
	//! The method to calculate current density
	void Calculate_J(const vector<vector<int> > &neighbour_list,
	           	     const vector<vector<double> > &coefficient_laplacian,
	                 const vector<vector<double> > &coefficient_dudx,
                     const vector<vector<double> > &coefficient_dudy,
                     const vector<vector<double> > &coefficient_dudz,
                     const vector<vector<double> > &coefficient_laplacian_boundary,
                     const vector<vector<double> > &coefficient_dudx_boundary,
                     const vector<vector<double> > &coefficient_dudy_boundary,
                     const vector<vector<double> > &coefficient_dudz_boundary,
                     const vector<double> &normal_x,
                     const vector<double> &normal_y,
                     const vector<double> &normal_z,
                     const vector<double> &Jx,
                     const vector<double> &Jy,
                     const vector<double> &Jz,
                     const vector<double> &Phi);

private:
	vector<double> m_UxB_x, m_UxB_y, m_UxB_z;
};

} /* namespace std */
#endif /* MHDSOLVER_H_ */
