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
#include "DATA.h"
#include "VerletScheme.h"
namespace std {

class Controller {
public:
	Controller();
	virtual ~Controller();
	//! The method to do initialization for controller
	void Initialization();
	//! Start the controller
	void Start();
	double Get_dt();
private:
	DATA m_data;
	Time_Integrator_Base *m_time_integrator;
	double m_time, m_dt, m_print_interval, m_max_time, m_next_print_time;
	int m_step, m_max_step, m_i_output, m_nprint;
};

} /* namespace std */
#endif /* CONTROLLER_H_ */
