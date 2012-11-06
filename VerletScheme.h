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
private:
	//! pointer to data
	DATA *m_data;
	Elliptic_Solver_Base *m_ES;
	Hyperbolic_Solver_Base *m_HS;
	double m_old_dt;
};

} /* namespace std */
#endif /* VERLETSCHEME_H_ */
