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
	Set_xp_current(m_data->Get_x());
	Set_yp_current(m_data->Get_y());
	Set_zp_current(m_data->Get_z());
	Set_up_current(m_data->Get_u());
	Set_vp_current(m_data->Get_v());
	Set_wp_current(m_data->Get_w());
	Set_rho_current(m_data->Get_rho());
	Set_energy_current(m_data->Get_energy());
	Set_pressure_current(m_data->Get_pressure());
}

void Verlet_Scheme::Integrate(const double &dt){
	m_ES->Solve();
	m_data->Print(0,0,NULL);
	assert(0);
	Set_old_current_states();
	m_HS->Set_dt(dt+m_old_dt);
	m_HS->Solve(
			m_data->Get_neighbour_list(),
			m_data->Get_coefficient_laplacian(),
			m_data->Get_coefficient_dudx(),
			m_data->Get_coefficient_dudy(),
			m_data->Get_coefficient_dudz(),
			m_data->Get_force_x(),
			m_data->Get_force_y(),
			m_data->Get_force_z(),
			m_xp_old,
			m_yp_old,
			m_zp_old,
			m_up_old,
			m_vp_old,
			m_wp_old,
			m_rho_old,
			m_energy_old,
			m_pressure_old,
			m_xp_current,
			m_yp_current,
			m_zp_current,
			m_up_current,
			m_vp_current,
			m_wp_current,
			m_rho_current,
			m_energy_current,
			m_pressure_current,
			//! output
			m_data->Get_x(),
			m_data->Get_y(),
			m_data->Get_z(),
			m_data->Get_u(),
			m_data->Get_v(),
			m_data->Get_w(),
			m_data->Get_rho(),
			m_data->Get_energy(),
			m_data->Get_pressure());
	Set_old_dt(dt);
}

void Verlet_Scheme::Set_old_current_states(){
	Set_xp_old(m_xp_current);
	Set_yp_old(m_yp_current);
	Set_zp_old(m_zp_current);
	Set_up_old(m_up_current);
	Set_vp_old(m_vp_current);
	Set_wp_old(m_wp_current);
	Set_rho_old(m_rho_current);
	Set_energy_old(m_energy_current);
	Set_pressure_old(m_pressure_current);
	Set_xp_current(m_data->Get_x());
	Set_yp_current(m_data->Get_y());
	Set_zp_current(m_data->Get_z());
	Set_up_current(m_data->Get_u());
	Set_vp_current(m_data->Get_v());
	Set_wp_current(m_data->Get_w());
	Set_rho_current(m_data->Get_rho());
	Set_energy_current(m_data->Get_energy());
	Set_pressure_current(m_data->Get_pressure());
}
} /* namespace std */
