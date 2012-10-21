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
namespace std {

class Verlet_Scheme: public std::Time_Integrator_Base {
public:
	Verlet_Scheme(DATA *data);
	virtual ~Verlet_Scheme();
	void Init();
	virtual void Integrate(double dt);
private:
	//! Variables of n - 1 step
	vector<double> m_xp_pre;
	vector<double> m_yp_pre;
	vector<double> m_zp_pre;
	vector<double> m_up_pre;
	vector<double> m_vp_pre;
	vector<double> m_wp_pre;
	vector<double> m_rho_pre;
	vector<double> m_energy_pre;
	//! pointer to data
	DATA *m_data;
	double m_dt;
	Elliptic_Solver_Base *m_ES;
	Hyperbolic_Solver_Base *m_HS;
};

} /* namespace std */
#endif /* VERLETSCHEME_H_ */
