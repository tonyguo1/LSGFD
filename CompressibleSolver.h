/*
 * CompressibleSolver.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The class to solve hyperbolic step for compressible material
 */
#ifndef COMPRESSIBLESOLVER_H_
#define COMPRESSIBLESOLVER_H_

#include "HyperbolicSolverBase.h"
#include "DATA.h"
#include "EOS.h"

namespace std {

class Compressible_Solver: public std::Hyperbolic_Solver_Base {
public:
	Compressible_Solver(DATA *data);
	virtual ~Compressible_Solver();
	//! The method to solve hyperbolic step
	virtual void Solve();
	virtual void Set_dt(double dt){m_dt = dt;};
	//! The method to update states
	/*!
	 *  This method will update velocity, density, energy, and position
	 *  For example: U_new = U_old + dt * F;
	 */
	void Update_States(//! input
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
			vector<double> &p_new);
	//! The methods to set everything
	void Set_xp_old(vector<double> &xp){m_xp_old.assign(xp.begin(),xp.end());}
	void Set_yp_old(vector<double> &yp){m_yp_old.assign(yp.begin(),yp.end());}
	void Set_zp_old(vector<double> &zp){m_zp_old.assign(zp.begin(),zp.end());}
	void Set_up_old(vector<double> &up){m_up_old.assign(up.begin(),up.end());}
	void Set_vp_old(vector<double> &vp){m_vp_old.assign(vp.begin(),vp.end());}
	void Set_wp_old(vector<double> &wp){m_wp_old.assign(wp.begin(),wp.end());}
	void Set_rho_old(vector<double> &rho){m_rho_old.assign(rho.begin(),rho.end());}
	void Set_energy_old(vector<double> &energy){m_energy_old.assign(energy.begin(),energy.end());}
	void Set_pressure_old(vector<double> &pressure){m_pressure_old.assign(pressure.begin(),pressure.end());}
	void Set_xp_current(vector<double> &xp){m_xp_current.assign(xp.begin(),xp.end());}
	void Set_yp_current(vector<double> &yp){m_yp_current.assign(yp.begin(),yp.end());}
	void Set_zp_current(vector<double> &zp){m_zp_current.assign(zp.begin(),zp.end());}
	void Set_up_current(vector<double> &up){m_up_current.assign(up.begin(),up.end());}
	void Set_vp_current(vector<double> &vp){m_vp_current.assign(vp.begin(),vp.end());}
	void Set_wp_current(vector<double> &wp){m_wp_current.assign(wp.begin(),wp.end());}
	void Set_rho_current(vector<double> &rho){m_rho_current.assign(rho.begin(),rho.end());}
	void Set_energy_current(vector<double> &energy){m_energy_current.assign(energy.begin(),energy.end());}
	void Set_pressure_current(vector<double> &pressure){m_pressure_current.assign(pressure.begin(),pressure.end());}
private:
	DATA *m_data;
	EOS_BASE *m_EOS;
	double m_dt;
	vector<double> m_xp_old;
	vector<double> m_yp_old;
	vector<double> m_zp_old;
	vector<double> m_up_old;
	vector<double> m_vp_old;
	vector<double> m_wp_old;
	vector<double> m_rho_old;
	vector<double> m_energy_old;
	vector<double> m_pressure_old;
	vector<double> m_xp_current;
	vector<double> m_yp_current;
	vector<double> m_zp_current;
	vector<double> m_up_current;
	vector<double> m_vp_current;
	vector<double> m_wp_current;
	vector<double> m_rho_current;
	vector<double> m_energy_current;
	vector<double> m_pressure_current;
};



} /* namespace std */
#endif /* COMPRESSIBLESOLVER_H_ */
