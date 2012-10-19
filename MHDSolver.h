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
	virtual void SetMatrix();
	//! The method to initialize UxB
	void Initializae_UxB();

private:
	vector<double> m_UxB_x, m_UxB_y, m_UxB_z;
};

} /* namespace std */
#endif /* MHDSOLVER_H_ */
