/*
 * MHDSolver.cpp
 *
 *  Created on: Oct 19, 2012
 *      Author: tongfei
 */

#include "MHDSolver.h"

namespace std {

MHD_Solver::MHD_Solver(DATA *data):Elliptic_Solver_Base(data) {
	// TODO Auto-generated constructor stub

}

MHD_Solver::~MHD_Solver() {
	// TODO Auto-generated destructor stub
}

virtual void MHD_Solver::Solve(){

}

virtual void MHD_Solver::SetMatrix(){

}

void MHD_Solver::Initializae_UxB(){
	vector<double> velo(3,0), UxB(3,0), B(3,0);
	m_UxB_x.assign(m_data->m_num_of_par,0);
	m_UxB_y.assign(m_data->m_num_of_par,0);
	m_UxB_z.assign(m_data->m_num_of_par,0);
	for (int i = 0; i < m_data->m_num_of_par; i++){
		m_data->Get_MagneticFiled(m_data->m_xp[i],m_data->m_yp[i],m_data->m_zp[i],B);
		velo[0] = m_data->m_up[i];
		velo[1] = m_data->m_vp[i];
		velo[2] = m_data->m_wp[i];
		m_data->Cross_Product(velo,B,UxB);
		m_UxB_x[i] = UxB[0];
		m_UxB_y[i] = UxB[1];
		m_UxB_z[i] = UxB[2];
	}
}
} /* namespace std */
