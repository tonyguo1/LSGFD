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

class Hyperbolic_Solver_Base {
public:
	//! The constructor, initial m_dt to be 0
	Hyperbolic_Solver_Base();
	virtual ~Hyperbolic_Solver_Base();
	//! The method to solve hyperbolic step
	virtual void Solve(){};
};

} /* namespace std */
#endif /* HYPERBOLICSOLVERBASE_H_ */
