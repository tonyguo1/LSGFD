/*
 * CompressibleSolver.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The class to solve hyperbolic step for compressible material
 */
#ifndef COMPRESSIBLESOLVER_H_
#define COMPRESSIBLESOLVER_H_

#include "HyperbolicSolverBase.h"
#include "DATA.h"

namespace std {

class Compressible_Solver: public std::Hyperbolic_Solver_Base {
public:
	Compressible_Solver(DATA *data);
	virtual ~Compressible_Solver();
	//! The method to solve hyperbolic step
	virtual void Set_dt(double dt){m_dt = dt;};
	//! The method to update states
	/*!
	 *  This method will update velocity, density, energy, and position
	 *  For example: U_new = U_old + dt * F;
	 */
	virtual void Solve(//! input
			const vector<vector<int> > &neighbour_list,
			const vector<vector<double> > &coefficient_laplacian,
			const vector<vector<double> > &coefficient_dudx,
			const vector<vector<double> > &coefficient_dudy,
			const vector<vector<double> > &coefficient_dudz,
			const vector<double> &force_x,
			const vector<double> &force_y,
			const vector<double> &force_z,
			const vector<double> &xp_old,
			const vector<double> &yp_old,
			const vector<double> &zp_old,
			const vector<double> &up_old,
			const vector<double> &vp_old,
			const vector<double> &wp_old,
			const vector<double> &rho_old,
			const vector<double> &e_old,
			const vector<double> &p_old,
			const vector<double> &xp_current,
			const vector<double> &yp_current,
			const vector<double> &zp_current,
			const vector<double> &up_current,
			const vector<double> &vp_current,
			const vector<double> &wp_current,
			const vector<double> &rho_current,
			const vector<double> &e_current,
			const vector<double> &p_current,
			//! output
			vector<double> &xp_new,
			vector<double> &yp_new,
			vector<double> &zp_new,
			vector<double> &up_new,
			vector<double> &vp_new,
			vector<double> &wp_new,
			vector<double> &rho_new,
			vector<double> &e_new,
			vector<double> &p_new);
private:
	DATA *m_data;
	double m_dt;
};



} /* namespace std */
#endif /* COMPRESSIBLESOLVER_H_ */
