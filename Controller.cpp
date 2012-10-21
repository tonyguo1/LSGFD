/*
 * Controller.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "Controller.h"

namespace std {

Controller::Controller():m_time(0),m_dt(0),m_step(0) {
	// TODO Auto-generated constructor stub
	m_time_integrator = NULL;

}

Controller::~Controller() {
	// TODO Auto-generated destructor stub
}

void Controller::Initialization(){
	m_data.Initialization();
    //! Initialization for time integrator
}

void Controller::Start(){
	Initialization();
	double max_time;
	int max_step;
	while (m_time < max_time && m_step < max_step)
	{
		//! Get dt according to states and print time
		m_dt = Get_dt();
		m_time += m_dt;
		m_step++;
		m_time_integrator->Integrate(m_dt);
		Print_control(m_time);
	}
	if (time < max_time)
		m_data.Print(m_time);
}

} /* namespace std */
