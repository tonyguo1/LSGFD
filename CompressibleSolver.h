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
};

} /* namespace std */
#endif /* COMPRESSIBLESOLVER_H_ */
