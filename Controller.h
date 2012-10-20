/*
 * Controller.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The class to do scheduling
 */

#ifndef CONTROLLER_H_
#define CONTROLLER_H_
#include <DATA.h>
namespace std {

class Controller {
public:
	Controller();
	virtual ~Controller();
	void start();
private:
	DATA *m_data;
	Elliptic_Solver_Base m_mhd_solver;
	Hyperbolic_Solver_Base m_compressible_solver;
};

} /* namespace std */
#endif /* CONTROLLER_H_ */
