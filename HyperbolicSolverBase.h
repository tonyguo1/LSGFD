/*
 * HyperbolicSolverbase.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The base class for all hyperbolic solver
 */
#ifndef HYPERBOLICSOLVERBASE_H_
#define HYPERBOLICSOLVERBASE_H_
#include <vector>
namespace std {

class Hyperbolic_Solver_Base {
public:
	//! The constructor, initial m_dt to be 0
	Hyperbolic_Solver_Base();
	virtual ~Hyperbolic_Solver_Base();
	//! The method to solve hyperbolic step
	virtual void Solve()=0;
	virtual void Set_dt(double dt)=0;
	//! The methods to set everything
	virtual void Set_xp_old(vector<double> &xp)=0;
	virtual void Set_yp_old(vector<double> &yp)=0;
	virtual void Set_zp_old(vector<double> &zp)=0;
	virtual void Set_up_old(vector<double> &up)=0;
	virtual void Set_vp_old(vector<double> &vp)=0;
	virtual void Set_wp_old(vector<double> &wp)=0;
	virtual void Set_rho_old(vector<double> &rho)=0;
	virtual void Set_energy_old(vector<double> &energy)=0;
	virtual void Set_pressure_old(vector<double> &pressure)=0;
	virtual void Set_xp_current(vector<double> &xp)=0;
	virtual void Set_yp_current(vector<double> &yp)=0;
	virtual void Set_zp_current(vector<double> &zp)=0;
	virtual void Set_up_current(vector<double> &up)=0;
	virtual void Set_vp_current(vector<double> &vp)=0;
	virtual void Set_wp_current(vector<double> &wp)=0;
	virtual void Set_rho_current(vector<double> &rho)=0;
	virtual void Set_energy_current(vector<double> &energy)=0;
	virtual void Set_pressure_current(vector<double> &pressure)=0;
	virtual vector<double>& Get_xp_current()=0;
	virtual vector<double>& Get_yp_current()=0;
	virtual vector<double>& Get_zp_current()=0;
	virtual vector<double>& Get_up_current()=0;
	virtual vector<double>& Get_vp_current()=0;
	virtual vector<double>& Get_wp_current()=0;
	virtual vector<double>& Get_rho_current()=0;
	virtual vector<double>& Get_energy_current()=0;
	virtual vector<double>& Get_pressure_current()=0;
};

} /* namespace std */
#endif /* HYPERBOLICSOLVERBASE_H_ */
