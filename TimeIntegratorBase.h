/*
 * TimeIntegratorBase.h
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */

#ifndef TIMEINTEGRATORBASE_H_
#define TIMEINTEGRATORBASE_H_

namespace std {

class Time_Integrator_Base {
public:
	//! Constructor with time step
	Time_Integrator_Base(double dt);
	virtual ~Time_Integrator_Base();
	virtual void Integrate() = 0;
private:
	double m_dt;
};

} /* namespace std */
#endif /* TIMEINTEGRATORBASE_H_ */
