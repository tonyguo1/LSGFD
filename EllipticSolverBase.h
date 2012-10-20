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
	//! The method to solve elliptic problem
	virtual void Solve()=0;
	//! The method for setting Matrix
	virtual void SetMatrix()=0;
private:
	//
	DATA *m_data;
	PETSc *m_petsc;
};

} /* namespace std */
#endif /* ELLIPTICSOLVERBASE_H_ */
