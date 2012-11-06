/*
 * CompressibleSolver.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "CompressibleSolver.h"
#include "assert.h"

namespace std {
Compressible_Solver::Compressible_Solver(DATA *data):m_data(data),m_dt(0){
}

Compressible_Solver::~Compressible_Solver() {
	// TODO Auto-generated destructor stub
}

void Compressible_Solver::Solve()
{
	Update_States(m_data->Get_neighbour_list(),
			m_data->Get_coefficient_laplacian(),
			m_data->Get_coefficient_dudx(),
			m_data->Get_coefficient_dudy(),
			m_data->Get_coefficient_dudz(),
			m_data->Get_force_x(),
			m_data->Get_force_y(),
			m_data->Get_force_z(),
			m_data->Get_x(),
			m_data->Get_y(),
			m_data->Get_z(),
			m_data->Get_u(),
			m_data->Get_v(),
			m_data->Get_w(),
			m_data->Get_rho(),
			m_data->Get_energy(),
			m_data->Get_pressure()
	);
}

void Compressible_Solver::Update_States(//! input
		const vector<vector<int> > &neighbour_list,
		const vector<vector<double> > &coefficient_laplacian,
		const vector<vector<double> > &coefficient_dudx,
		const vector<vector<double> > &coefficient_dudy,
		const vector<vector<double> > &coefficient_dudz,
		const vector<double> &force_x,
		const vector<double> &force_y,
		const vector<double> &force_z,
		//! output
		vector<double> &xp_new,
		vector<double> &yp_new,
		vector<double> &zp_new,
		vector<double> &up_new,
		vector<double> &vp_new,
		vector<double> &wp_new,
		vector<double> &rho_new,
		vector<double> &e_new,
		vector<double> &p_new){
	double cmax = 0;
	int num_of_par = m_data->Get_num_of_par();
	double ax,ay,az,arho,ae;
	for (int i_index = 0; i_index < num_of_par; i_index++){
		double cs = 0;
		int num_neigh = neighbour_list[i_index].size();
		ax = ay = az = arho = ae = 0;
		for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
			int nei_index = neighbour_list[i_index][i_neigh];
			ax   += m_pressure_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			ay   += m_pressure_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			az   += m_pressure_current[nei_index]*coefficient_dudz[i_index][i_neigh];
			arho += m_up_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			arho += m_vp_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			arho += m_wp_current[nei_index]*coefficient_dudz[i_index][i_neigh];
			ae   += m_up_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			ae   += m_vp_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			ae   += m_wp_current[nei_index]*coefficient_dudz[i_index][i_neigh];
		}
		if (m_rho_current[i_index]!=0){
			ax   = -ax/m_rho_current[i_index] + force_x[i_index];
			ay   = -ay/m_rho_current[i_index] + force_y[i_index];
			az   = -az/m_rho_current[i_index] + force_z[i_index];
			ae   = -ae*m_pressure_current[i_index]/m_rho_current[i_index];
			arho = -arho*m_rho_current[i_index];
		}
		else {assert(0);}

		xp_new[i_index]  = m_xp_old[i_index]     + m_up_current[i_index]*m_dt;
		yp_new[i_index]  = m_yp_old[i_index]     + m_vp_current[i_index]*m_dt;
		zp_new[i_index]  = m_zp_old[i_index]     + m_wp_current[i_index]*m_dt;
		up_new[i_index]  = m_up_old[i_index]     + ax*m_dt;
		vp_new[i_index]  = m_vp_old[i_index]     + ay*m_dt;
		wp_new[i_index]  = m_zp_old[i_index]     + az*m_dt;
		rho_new[i_index] = m_rho_old[i_index]    + arho*m_dt;
		e_new[i_index]   = m_energy_old[i_index] + ae*m_dt;
		m_data->Set_eos(rho_new[i_index], e_new[i_index], p_new[i_index], cs);
		if (cs > cmax)
			cmax = cs;
	}
	m_data->Set_cmax(cmax);
}

} /* namespace std */
