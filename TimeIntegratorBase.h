/*
 * TimeIntegratorBase.h
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */

#ifndef TIMEINTEGRATORBASE_H_
#define TIMEINTEGRATORBASE_H_
#include "Solver.h"
/**
 *
 * \brief	The base class of all time integrator
 */
namespace std {

class Time_Integrator_Base {
public:
	//! Constructor with time step
	Time_Integrator_Base();
	virtual ~Time_Integrator_Base();
	virtual void Integrate(double dt) = 0;
};

} /* namespace std */
#endif /* TIMEINTEGRATORBASE_H_ */
