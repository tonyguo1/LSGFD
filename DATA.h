/*
 * DATA.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
//! The class to store basic data shared by all solvers and some data manipulation method
/*!
 *  The pointer to this class will be pass to solvers by controller
 */
#ifndef DATA_H_
#define DATA_H_
#include <vector>
#include <particles.h>
#include <deque>
#include <utility>
#include <algorithm>
namespace std {

class DATA {
public:
	DATA();
	virtual ~DATA();
	//! The member function to do initialization
	/*!
	 *  Initialize particles and some settings
	 */
	void Initialization();
	//! The member function to build neighbour list and coefficient list
	/*!
	 *  This function will call functions from octree search class and coefficient calculation class
	 */
	void Buildup_neigh_list_and_ceoff();
	//! The member function to clear neighbour list and coefficient list
	void clear_neigh_list_and_ceoff();


private:
	//! Particle number
	size_t m_num_of_par;
	//! initial shortest particle distance
	double m_distance;
    //! Particle position and velocity, density, volume, mass
    vector<double> m_xp, m_yp, m_zp, m_up, m_vp, m_wp, m_rho, m_vol, m_mass;
    //! Compressible code related quantities, energy density, pressure, temperature
    vector<double> m_e, m_p, m_T;
    //! MHD related quantities
    vecotr<double> m_Jx, m_Jy, m_Jz, m_phi;
    //! The neighbour list
    /*!
     *  A vector of vector contains all the neighbours of all particles, the first place is the particle itself
     *  should be set in each time step
     */
    vector<vector<int> > m_neighbour_list;
    //! The coefficient list
    /*!
     *  Store the coefficients of dudx, dudy, dudz, and laplacian of all particles, each coefficient is for the particle in the corresponding place of neighbour_list
     */
    vector<vector<double> > m_coefficient_laplacian, m_coefficient_dudx, m_coefficient_dudy, m_coefficient_dudz;
};

} /* namespace std */
#endif /* DATA_H_ */
