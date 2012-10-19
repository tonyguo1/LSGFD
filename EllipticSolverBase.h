/*
 * EllipticSolverBase.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
//! The base class for all elliptic solver
/*!
 *  This class will be inherited by all elliptic solvers.
 */

#ifndef ELLIPTICSOLVERBASE_H_
#define ELLIPTICSOLVERBASE_H_
#include <DATA.h>
#include <solver.h>

namespace std {

class Elliptic_Solver_Base {
public:
	//! Constructor without any parameter
	/*!
	 *  The basic one
	 */
	Elliptic_Solver_Base(DATA *data);
	//! Basic destructor
	virtual ~Elliptic_Solver_Base();
	//! The member function to solve elliptic problem
	virtual void Solve(){};
	//! The member function for setting Matrix
	/*!
	 *  For base form, no parameter and no return value
	 */
	virtual void SetMatrix(){};
	//! The member funtion for solving matrix
	/*!
	 *  For base form, no parameter and
	 */
private:
	//
	DATA *m_data;
	PETSc *m_petsc;
};

} /* namespace std */
#endif /* ELLIPTICSOLVERBASE_H_ */
