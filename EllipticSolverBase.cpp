/*
 * EllipticSolverBase.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "EllipticSolverBase.h"

namespace std {

Elliptic_Solver_Base::Elliptic_Solver_Base(DATA *data):m_data(data) {
	// TODO Auto-generated constructor stub
	m_petsc = NULL;
}

Elliptic_Solver_Base::~Elliptic_Solver_Base() {
	// TODO Auto-generated destructor stub
	if (m_petsc != NULL)
		delete m_petsc;
}

} /* namespace std */
