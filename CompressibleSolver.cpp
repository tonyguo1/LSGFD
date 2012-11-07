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

void Compressible_Solver::Solve(//! input
		const vector<vector<int> > &neighbour_list,
		const vector<vector<double> > &coefficient_laplacian,
		const vector<vector<double> > &coefficient_dudx,
		const vector<vector<double> > &coefficient_dudy,
		const vector<vector<double> > &coefficient_dudz,
		const vector<double> &force_x,
		const vector<double> &force_y,
		const vector<double> &force_z,
		const vector<double> &xp_old,
		const vector<double> &yp_old,
		const vector<double> &zp_old,
		const vector<double> &up_old,
		const vector<double> &vp_old,
		const vector<double> &wp_old,
		const vector<double> &rho_old,
		const vector<double> &energy_old,
		const vector<double> &pressure_old,
		const vector<double> &xp_current,
		const vector<double> &yp_current,
		const vector<double> &zp_current,
		const vector<double> &up_current,
		const vector<double> &vp_current,
		const vector<double> &wp_current,
		const vector<double> &rho_current,
		const vector<double> &energy_current,
		const vector<double> &pressure_current,
		//! output
		vector<double> &xp_new,
		vector<double> &yp_new,
		vector<double> &zp_new,
		vector<double> &up_new,
		vector<double> &vp_new,
		vector<double> &wp_new,
		vector<double> &rho_new,
		vector<double> &energy_new,
		vector<double> &pressure_new){
	Set_force(m_data->Get_x(),m_data->Get_y(),m_data->Get_z(),m_data->Get_force_x(),m_data->Get_force_y(),m_data->Get_force_z());
	double cmax = 0;
	int num_of_par = m_data->Get_num_of_par();
	double ax,ay,az,arho,ae;
	for (int i_index = 0; i_index < num_of_par; i_index++){
		double cs = 0;
		int num_neigh = neighbour_list[i_index].size();
		ax = ay = az = arho = ae = 0;
		for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
			int nei_index = neighbour_list[i_index][i_neigh];
			ax   += pressure_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			ay   += pressure_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			az   += pressure_current[nei_index]*coefficient_dudz[i_index][i_neigh];
			arho += up_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			arho += vp_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			arho += wp_current[nei_index]*coefficient_dudz[i_index][i_neigh];
			ae   += up_current[nei_index]*coefficient_dudx[i_index][i_neigh];
			ae   += vp_current[nei_index]*coefficient_dudy[i_index][i_neigh];
			ae   += wp_current[nei_index]*coefficient_dudz[i_index][i_neigh];
		}
		if (rho_current[i_index]!=0){
			ax   = -ax/rho_current[i_index] + force_x[i_index];
			ay   = -ay/rho_current[i_index] + force_y[i_index];
			az   = -az/rho_current[i_index] + force_z[i_index];
			ae   = -ae*pressure_current[i_index]/rho_current[i_index];
			arho = -arho*rho_current[i_index];
		}
		else {assert(0);}

		xp_new[i_index]  = xp_old[i_index]     + up_current[i_index]*m_dt;
		yp_new[i_index]  = yp_old[i_index]     + vp_current[i_index]*m_dt;
		zp_new[i_index]  = zp_old[i_index]     + wp_current[i_index]*m_dt;
		up_new[i_index]  = up_old[i_index]     + ax*m_dt;
		vp_new[i_index]  = vp_old[i_index]     + ay*m_dt;
		wp_new[i_index]  = wp_old[i_index]     + az*m_dt;
		rho_new[i_index] = rho_old[i_index]    + arho*m_dt;
		energy_new[i_index]   = energy_old[i_index] + ae*m_dt;
		m_data->Set_eos(rho_new[i_index], energy_new[i_index], pressure_new[i_index], cs);
		if (cs > cmax)
			cmax = cs;
	}
	m_data->Set_cmax(cmax);
}

void Compressible_Solver::Set_force(const vector<double> &x, const vector<double> &y, const vector<double> &z, vector<double> &fx, vector<double> &fy, vector<double> &fz){
	int num_of_par = m_data->Get_num_of_par();
	for (int i_index = 0; i_index < num_of_par; i_index++){
		if (z[i_index] > 0){
			fx[i_index] = 0;
			fy[i_index] = 0;
			fz[i_index] = 0;
		}
		else{
			fx[i_index] = 0;
			fy[i_index] = 0;
			fz[i_index] = 0;
		}
	}
}

} /* namespace std */
