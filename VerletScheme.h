/*
 * VerletScheme.h
 *
 *  Created on: Oct 20, 2012
 *      Author: tongfei
 */

#ifndef VERLETSCHEME_H_
#define VERLETSCHEME_H_

#include "TimeIntegratorBase.h"
#include <vector>
namespace std {

class Verlet_Scheme: public std::Time_Integrator_Base {
public:
	Verlet_Scheme();
	virtual ~Verlet_Scheme();
	virtual void Integrate();
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
};

} /* namespace std */
#endif /* VERLETSCHEME_H_ */
