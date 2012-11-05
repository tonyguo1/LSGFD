/*
 * MHDSolver.cpp
 *
 *  Created on: Oct 19, 2012
 *      Author: tongfei
 */

#include "MHDSolver.h"
#include <assert.h>
#include <iostream>
#include <stdio.h>
namespace std {

MHD_Solver::MHD_Solver(DATA *data):m_data(data) {
	// TODO Auto-generated constructor stub
	m_petsc = NULL;
	m_data->MHD_init();
}

MHD_Solver::~MHD_Solver() {
	// TODO Auto-generated destructor stub
}

void MHD_Solver::Solve(){
	PETSc petsc;
	m_petsc = &petsc;
	int num_of_par = m_data->Get_num_of_par();
	m_petsc->Create(0, num_of_par - 1, 30, 0);
	m_petsc->SetTol(1e-15);
	Initializae_UxB();
	SetMatrix(m_data->Get_neighbour_list(),
			  m_data->Get_coefficient_laplacian(),
			  m_data->Get_coefficient_dudx(),
			  m_data->Get_coefficient_dudy(),
			  m_data->Get_coefficient_dudz(),
			  m_data->Get_coefficient_laplacian_boundary(),
			  m_data->Get_coefficient_dudx_boundary(),
			  m_data->Get_coefficient_dudy_boundary(),
			  m_data->Get_coefficient_dudz_boundary(),
			  m_data->Get_normal_x(),
			  m_data->Get_normal_y(),
			  m_data->Get_normal_z());
	//m_petsc->Print_A(NULL);
	//m_petsc->Print_b(NULL);
	//m_petsc->Solve_withPureNeumann_GMRES();
	m_petsc->Solve_LSQR();
	//m_petsc->Print_x(NULL);
	int iter;
	double residual;
	m_petsc->GetNumIterations(&iter);
	m_petsc->GetFinalRelativeResidualNorm(&residual);
	double *x = new double[num_of_par];
	m_petsc->Get_x(x);
	Set_Phi(m_data->Get_Phi(), x);
	delete [] x;
	Calculate_J(m_data->Get_neighbour_list(),
		      	m_data->Get_coefficient_laplacian(),
		      	m_data->Get_coefficient_dudx(),
		      	m_data->Get_coefficient_dudy(),
		      	m_data->Get_coefficient_dudz(),
		      	m_data->Get_coefficient_laplacian_boundary(),
		      	m_data->Get_coefficient_dudx_boundary(),
		      	m_data->Get_coefficient_dudy_boundary(),
		      	m_data->Get_coefficient_dudz_boundary(),
		        m_data->Get_normal_x(),
		     	m_data->Get_normal_y(),
		      	m_data->Get_normal_z(),
		      	m_data->Get_Jx(),
		      	m_data->Get_Jy(),
		      	m_data->Get_Jz(),
		      	m_data->Get_Phi());
}

void MHD_Solver::SetMatrix(const vector<vector<int> > &neighbour_list,
		                           const vector<vector<double> > &coefficient_laplacian,
		                           const vector<vector<double> > &coefficient_dudx,
		                           const vector<vector<double> > &coefficient_dudy,
		                           const vector<vector<double> > &coefficient_dudz,
		                           const vector<vector<double> > &coefficient_laplacian_boundary,
		                           const vector<vector<double> > &coefficient_dudx_boundary,
		                           const vector<vector<double> > &coefficient_dudy_boundary,
		                           const vector<vector<double> > &coefficient_dudz_boundary,
		                           const vector<double> &normal_x,
		                           const vector<double> &normal_y,
		                           const vector<double> &normal_z){
	int boundary_index = 0;
	int num_of_par = m_data->Get_num_of_par();
	for (int index = 0; index < num_of_par; index++){
		int num_neigh = neighbour_list[index].size();
		double rhs = 0;
		if (m_data->Get_Boundary_Flag(index) == 1){
			for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
				double value = coefficient_laplacian_boundary[boundary_index][i_neigh];
				value = value - coefficient_laplacian_boundary[boundary_index][num_neigh] *
				(normal_x[index] * coefficient_dudx_boundary[boundary_index][i_neigh] + normal_y[index] * coefficient_dudy_boundary[boundary_index][i_neigh] + normal_z[index] * coefficient_dudz_boundary[boundary_index][i_neigh]) /
				(normal_x[index] * coefficient_dudx_boundary[boundary_index][num_neigh] + normal_y[index] * coefficient_dudy_boundary[boundary_index][num_neigh] + normal_z[index] * coefficient_dudz_boundary[boundary_index][num_neigh]);
				m_petsc->Add_A(index, neighbour_list[index][i_neigh], value);
				rhs += coefficient_dudx[index][i_neigh] * m_UxB_x[index] + coefficient_dudy[index][i_neigh] * m_UxB_y[index] + coefficient_dudz[index][i_neigh] * m_UxB_z[index];
			}
			double value =  normal_x[index] * m_UxB_x[index] + normal_y[index] * m_UxB_y[index] + normal_z[index] * m_UxB_z[index];
			rhs -= value * coefficient_laplacian_boundary[boundary_index][num_neigh] /
					(normal_x[index] * coefficient_dudx_boundary[boundary_index][num_neigh] + normal_y[index] * coefficient_dudy_boundary[boundary_index][num_neigh] + normal_z[index] * coefficient_dudz_boundary[boundary_index][num_neigh]);
			m_petsc->Add_b(index, rhs);
			++boundary_index;;
		}
		else{
			int num_neigh = neighbour_list[index].size();
			double rhs = 0;
			for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
				double value = coefficient_laplacian[index][i_neigh];
				m_petsc->Add_A(index, neighbour_list[index][i_neigh], value);
				rhs += coefficient_dudx[index][i_neigh] * m_UxB_x[index] + coefficient_dudy[index][i_neigh] * m_UxB_y[index] + coefficient_dudz[index][i_neigh] * m_UxB_z[index];
			}
			m_petsc->Add_b(index, rhs);
		}
	}
}

void MHD_Solver::Initializae_UxB(){
	double lightSpeed=3e7;
	vector<double> velo(3,0), UxB(3,0), B(3,0);
	int num_of_par = m_data->Get_num_of_par();
	m_UxB_x.assign(num_of_par,0);
	m_UxB_y.assign(num_of_par,0);
	m_UxB_z.assign(num_of_par,0);
	for (int i = 0; i < num_of_par; i++){
		m_data->Get_MagneticFiled(m_data->Get_x(i),m_data->Get_y(i),m_data->Get_z(i),B);
		velo[0] = m_data->Get_u(i);
		velo[1] = m_data->Get_v(i);
		velo[2] = m_data->Get_w(i);
		m_data->Cross_Product(velo,B,UxB);
		m_UxB_x[i] = UxB[0] / lightSpeed;
		m_UxB_y[i] = UxB[1] / lightSpeed;
		m_UxB_z[i] = UxB[2] / lightSpeed;
	}
}

void MHD_Solver::Calculate_J(const vector<vector<int> > &neighbour_list,
                             const vector<vector<double> > &coefficient_laplacian,
                             const vector<vector<double> > &coefficient_dudx,
                             const vector<vector<double> > &coefficient_dudy,
                             const vector<vector<double> > &coefficient_dudz,
                             const vector<vector<double> > &coefficient_laplacian_boundary,
                             const vector<vector<double> > &coefficient_dudx_boundary,
                             const vector<vector<double> > &coefficient_dudy_boundary,
                             const vector<vector<double> > &coefficient_dudz_boundary,
                             const vector<double> &normal_x,
                             const vector<double> &normal_y,
                             const vector<double> &normal_z,
                             vector<double> &Jx,
                             vector<double> &Jy,
                             vector<double> &Jz,
                             const vector<double> &Phi){
	double m_fluidConductivity = 1e13;
	int boundary_index = 0;
	int num_of_par = m_data->Get_num_of_par();
	for (int index = 0; index < num_of_par; index++){
		int num_neigh = neighbour_list[index].size();
		if (m_data->Get_Boundary_Flag(index) == 1){
			double value = normal_x[index] * m_UxB_x[index] + normal_y[index] * m_UxB_y[index] + normal_z[index] * m_UxB_z[index];
			double ghost = value / (normal_x[index] * coefficient_dudx_boundary[boundary_index][num_neigh]
			                      + normal_y[index] * coefficient_dudy_boundary[boundary_index][num_neigh]
			                      + normal_z[index] * coefficient_dudz_boundary[boundary_index][num_neigh]);
			for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
				int index1 = neighbour_list[index][i_neigh];
				Jx[index] += coefficient_dudx_boundary[boundary_index][i_neigh] * Phi[index1];
				Jy[index] += coefficient_dudy_boundary[boundary_index][i_neigh] * Phi[index1];
				Jz[index] += coefficient_dudz_boundary[boundary_index][i_neigh] * Phi[index1];
				ghost -= (normal_x[index] * coefficient_dudx_boundary[boundary_index][i_neigh] + normal_y[index] * coefficient_dudy_boundary[boundary_index][i_neigh] + normal_z[index] * coefficient_dudz_boundary[boundary_index][i_neigh]) * Phi[index1] /
						 (normal_x[index] * coefficient_dudx_boundary[boundary_index][num_neigh] + normal_y[index] * coefficient_dudy_boundary[boundary_index][num_neigh] + normal_z[index] * coefficient_dudz_boundary[boundary_index][num_neigh]);
			}
			Jx[index] += coefficient_dudx_boundary[boundary_index][num_neigh] * ghost;
			Jy[index] += coefficient_dudy_boundary[boundary_index][num_neigh] * ghost;
			Jz[index] += coefficient_dudz_boundary[boundary_index][num_neigh] * ghost;

			Jx[index]= m_fluidConductivity * (-Jx[index] + m_UxB_x[index]);
			Jy[index]= m_fluidConductivity * (-Jy[index] + m_UxB_y[index]);
			Jz[index]= m_fluidConductivity * (-Jz[index] + m_UxB_z[index]);
			boundary_index++;
		}
		else{
			for (int i_neigh = 0; i_neigh < num_neigh; i_neigh++){
				int index1 = neighbour_list[index][i_neigh];
				Jx[index] += coefficient_dudx[index][i_neigh] * Phi[index1];
				Jy[index] += coefficient_dudy[index][i_neigh] * Phi[index1];
				Jz[index] += coefficient_dudz[index][i_neigh] * Phi[index1];
			}
			Jx[index]= m_fluidConductivity * (-Jx[index] + m_UxB_x[index]);
			Jy[index]= m_fluidConductivity * (-Jy[index] + m_UxB_y[index]);
			Jz[index]= m_fluidConductivity * (-Jz[index] + m_UxB_z[index]);
		}
	}
}

void MHD_Solver::Set_Phi(vector<double> &phi, double *x){
	int num_of_par = m_data->Get_num_of_par();
	for (int index = 0; index < num_of_par; index++)
		phi[index] = x[index];
}

void MHD_Solver::Set_force(const vector<double> &Jx, const vector<double> &Jy, const vector<double> &Jz, vector<double> &fx, vector<double> &fy, vector<double> &fz){
	int num_of_par = m_data->Get_num_of_par();
	double lightSpeed=3e7;
	for (int i = 0; i < num_of_par; i++){
		vector<double> J(3,0), B(3,0);
		vector<double> rslt(3,0);
		J[0] = Jx[i];
		J[1] = Jy[i];
		J[2] = Jz[i];
		m_data->Get_MagneticFiled(m_data->Get_x(i),m_data->Get_y(i),m_data->Get_z(i),B);
		m_data->Cross_Product(J,B,rslt);
		fx[i] = rslt[0] / lightSpeed;
		fy[i] = rslt[1] / lightSpeed;
		fz[i] = rslt[2] / lightSpeed;
	}
}
} /* namespace std */
