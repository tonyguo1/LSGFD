/*
 * VerletScheme.cpp
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */

#include "VerletScheme.h"
#include "assert.h"
namespace std {

Verlet_Scheme::Verlet_Scheme(DATA *data):m_data(data),m_old_dt(0) {
	Init();
}

Verlet_Scheme::~Verlet_Scheme() {
	// TODO Auto-generated destructor stub
}

void Verlet_Scheme::Init(){
    m_ES = new MHD_Solver(m_data);
    m_HS = new Compressible_Solver(m_data);
	m_HS->Set_xp_current(m_data->Get_x());
	m_HS->Set_yp_current(m_data->Get_y());
	m_HS->Set_zp_current(m_data->Get_z());
	m_HS->Set_up_current(m_data->Get_u());
	m_HS->Set_vp_current(m_data->Get_v());
	m_HS->Set_wp_current(m_data->Get_w());
	m_HS->Set_rho_current(m_data->Get_rho());
	m_HS->Set_energy_current(m_data->Get_energy());
	m_HS->Set_pressure_current(m_data->Get_pressure());
}

void Verlet_Scheme::Integrate(const double &dt){
	m_ES->Solve();
	m_data->Print(0,0,NULL);
	assert(0);
	m_HS->Set_xp_old(m_HS->Get_xp_current());
	m_HS->Set_yp_old(m_HS->Get_yp_current());
	m_HS->Set_zp_old(m_HS->Get_zp_current());
	m_HS->Set_up_old(m_HS->Get_up_current());
	m_HS->Set_vp_old(m_HS->Get_vp_current());
	m_HS->Set_wp_old(m_HS->Get_wp_current());
	m_HS->Set_rho_old(m_HS->Get_rho_current());
	m_HS->Set_energy_old(m_HS->Get_energy_current());
	m_HS->Set_pressure_old(m_HS->Get_pressure_current());
	m_HS->Set_xp_current(m_data->Get_x());
	m_HS->Set_yp_current(m_data->Get_y());
	m_HS->Set_zp_current(m_data->Get_z());
	m_HS->Set_up_current(m_data->Get_u());
	m_HS->Set_vp_current(m_data->Get_v());
	m_HS->Set_wp_current(m_data->Get_w());
	m_HS->Set_rho_current(m_data->Get_rho());
	m_HS->Set_energy_current(m_data->Get_energy());
	m_HS->Set_pressure_current(m_data->Get_pressure());
	m_HS->Set_dt(dt+m_old_dt);
	m_HS->Solve();
	m_old_dt = dt;
}

} /* namespace std */
