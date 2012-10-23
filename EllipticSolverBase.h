/*
 * EllipticSolverBase.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The base class for all elliptic solver
 */

/*!
 *  This class will be inherited by all elliptic solvers.
 */

#ifndef ELLIPTICSOLVERBASE_H_
#define ELLIPTICSOLVERBASE_H_
#include "DATA.h"
#include "solver_petsc.h"

namespace std {

class Elliptic_Solver_Base {
public:
	//! Constructor without any parameter
	/*!
	 *  The basic one
	 */
	Elliptic_Solver_Base();
	//! Basic destructor
	virtual ~Elliptic_Solver_Base();
	//! The method to solve elliptic problem
	virtual void Solve(){};
	//! The method for setting Matrix
};

} /* namespace std */
#endif /* ELLIPTICSOLVERBASE_H_ */
