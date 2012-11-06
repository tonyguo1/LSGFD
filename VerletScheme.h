/*
 * VerletScheme.h
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The class to perform Verlet Scheme
 */
#ifndef VERLETSCHEME_H_
#define VERLETSCHEME_H_

#include "TimeIntegratorBase.h"
#include <vector>
#include <Solver.h>
namespace std {

class Verlet_Scheme: public std::Time_Integrator_Base {
public:
	Verlet_Scheme(DATA *data);
	virtual ~Verlet_Scheme();
	void Init();
	virtual void Integrate(const double &dt);
	void Set_old_current_states();
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
	void Set_old_dt(double dt){m_old_dt = dt;}
private:
	//! pointer to data
	DATA *m_data;
	Elliptic_Solver_Base *m_ES;
	Hyperbolic_Solver_Base *m_HS;
	double m_old_dt;
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
#endif /* VERLETSCHEME_H_ */
