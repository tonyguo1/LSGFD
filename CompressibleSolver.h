/*
 * CompressibleSolver.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#ifndef COMPRESSIBLESOLVER_H_
#define COMPRESSIBLESOLVER_H_

#include "HyperbolicSolverbase.h"

namespace std {

class Compressible_Solver: public std::Hyperbolic_Solver_base {
public:
	Compressible_Solver();
	virtual ~Compressible_Solver();
	//! The method to solve hyperbolic step
	virtual void Solve();
	//! The method to update states
	/*!
	 *  This method will update velocity, density, energy, and position
	 *  For example: U_new = U_old + dt * F;
	 */
	void Update_States(//! input
				       const vector<vector<int> > &neighbour_list,
				       const vector<vector<double> > &coefficient_laplacian,
	           	       const vector<vector<double> > &coefficient_dudx,
	           	       const vector<vector<double> > &coefficient_dudy,
	           	       const vector<vector<double> > &coefficient_dudz,
	           	       const vector<double> &xp_old,
	           	       const vector<double> &yp_old,
	           	       const vector<double> &zp_old,
	           	       const vector<double> &up_old,
	           	       const vector<double> &vp_old,
	           	       const vector<double> &wp_old,
	           	       const vector<double> &rho_old,
	           	       const vector<double> &energy_old,
	           	       const vector<double> &pressure_mid,
	           	       const vector<double> &up_mid,
	           	       const vector<double> &vp_mid,
	           	       const vector<double> &wp_mid,
	           	       //! output
	           	       vector<double> &xp_new,
	           	       vector<double> &yp_new,
	           	       vector<double> &zp_new,
	           	       vector<double> &up_new,
	           	       vector<double> &vp_new,
	           	       vector<double> &wp_new,
	           	       vector<double> &rho_new,
	           	       vector<double> &e_new);
};

} /* namespace std */
#endif /* COMPRESSIBLESOLVER_H_ */
