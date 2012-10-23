/*
 * VerletScheme.cpp
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */

#include "VerletScheme.h"

namespace std {

Verlet_Scheme::Verlet_Scheme(DATA *data):m_data(data),m_dt(0.005) {
	Init();
}

Verlet_Scheme::~Verlet_Scheme() {
	// TODO Auto-generated destructor stub
}

void Verlet_Scheme::Init(){
	m_xp_pre = m_data->Get_x();
	m_yp_pre = m_data->Get_y();
	m_zp_pre = m_data->Get_z();
	m_up_pre = m_data->Get_u();
	m_vp_pre = m_data->Get_v();
	m_wp_pre = m_data->Get_w();
	m_rho_pre = m_data->Get_rho();
    m_energy_pre = m_data->Get_energy();
    m_ES = new MHD_Solver(m_data);
    m_HS = new Compressible_Solver(m_data);
}

void Verlet_Scheme::Integrate(const double &dt){
	m_ES->Solve();
	m_HS->Set_dt(dt);
	m_HS->Solve();
}

} /* namespace std */
