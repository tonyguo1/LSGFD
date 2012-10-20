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

namespace std {

class Hyperbolic_Solver_base {
public:
	//! The constructor, initial m_dt to be 0
	Hyperbolic_Solver_base();
	virtual ~Hyperbolic_Solver_base();
	//! The method to solve hyperbolic step
	virtual void Solve()=0;
	//! The method to set m_dt
	void set_dt(double dt){m_dt = dt;}
private:
	double m_dt;
};

} /* namespace std */
#endif /* HYPERBOLICSOLVERBASE_H_ */
