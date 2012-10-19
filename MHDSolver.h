/*
 * MHDSolver.h
 *
 *  Created on: Oct 19, 2012
 *      Author: tongfei
 */

#ifndef MHDSOLVER_H_
#define MHDSOLVER_H_

#include "EllipticSolverBase.h"

namespace std {

class MHD_Solver: public std::Elliptic_Solver_Base {
public:
	MHD_Solver();
	virtual ~MHD_Solver();
};

} /* namespace std */
#endif /* MHDSOLVER_H_ */
