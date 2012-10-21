/*
 * Controller.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "Controller.h"

namespace std {

Controller::Controller() {
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
	double max_time, time = 0, dt;
	int max_step, step;
	while (time < max_time && step < max_step)
	{
		//! Get dt according to states and print time
		dt = m_data.Get_dt();
		time += dt;
		step++;
		m_time_integrator->Integrate(dt);
		m_data.Print_control(time);
	}
	if (time < max_time)
		m_data.Print(time);
}

} /* namespace std */
