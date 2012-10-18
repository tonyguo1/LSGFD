/*
 * MHDsolver.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
//! The class to solve electric potential
/*!
 *  The governing equation is $\nabla ^2 \varphi = \nabla \cdot \frac{1c} (\mathbf{U} \times \mathbf{B})$
 *  The boundary condition is $\frac{\partial\varphi}{\partial n} = (\mathbf{U} \times \mathbf{B})\cdot\mathbf{n}$
 */
#ifndef MHDSOLVER_H_
#define MHDSOLVER_H_

#include "EllipticSolverBase.h"

namespace std {

class MHD_solver: public std::Elliptic_Solver_Base {
public:
	MHD_solver();
	virtual ~MHD_solver();
	//! The method to solve electric potential and calculate current density
	virtual void solve();

    //! The ethod to setup matrix
    /*!
     *  \param neighbour_list neighbour list of particles
     *  \param laplacian coefficients for laplacian operator
     *  \param coefficient_dudx coefficients for first derivative in x direction
     *  \param coefficient_dudy coefficients for first derivative in y direction
     *  \param coefficient_dudz coefficients for first derivative in z direction
     */
	virtual void SetMatrix(
			const vector<vector<int> > &neighbour_list,
			const vector<vector<int> > &coefficient_laplacian,
			const vector<vector<int> > &coefficient_dudx,
			const vector<vector<int> > &coefficient_dudy,
			const vector<vector<int> > &coefficient_dudz,
				);
};

} /* namespace std */
#endif /* MHDSOLVER_H_ */
