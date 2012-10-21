/*
 * Controller.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "Controller.h"

namespace std {

Controller::Controller():m_time(0),m_dt(0),m_step(0),m_i_output(0),m_nprint(0),m_max_step(0),m_max_time(0),m_next_print_time(0){
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
	while (m_time < m_max_time && m_step < m_max_step)
	{
		//! Get dt according to states and print time
		m_dt = Get_dt();
		m_time += m_dt;
		m_step++;
		m_time_integrator->Integrate(m_dt);
		if (m_i_output)
			m_data.Print(m_time);
	}
	if (time < m_max_time)
		m_data.Print(m_time);
}

double Controller::Get_dt(){
	double dt = m_data.Get_dt();
	while (m_step == 1 && m_next_print_time - m_time < dt)
		m_next_print_time += m_print_interval;
	if (m_time + dt >= m_next_print_time && m_time + dt < m_max_time)
	{
		dt = m_next_print_time - m_time + 1e-16;
		m_next_print_time += m_print_interval;
		if (m_next_print_time > m_max_time)
			m_next_print_time = m_max_time;
		m_i_output = 1;
		m_nprint++;
	}
	else if (m_time +dt >= m_max_time)
	{
		dt = m_max_time - m_time + 1e-16;
		m_i_output = 1;
		m_nprint++;
	}
	return dt;
}

} /* namespace std */
