/*
 * DATA.h
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */
/**
 *
 * \brief	The class to store basic data shared by all solvers and some data manipulation method
 */

/*!
 *  The pointer to this class will be pass to solvers by controller
 */
#ifndef DATA_H_
#define DATA_H_
#include <vector>
#include "octree.h"
#include "solver_lapack_cf.h"
#include <deque>
#include <utility>
#include <algorithm>
#define PI 3.14159265359
namespace std {

class DATA {
public:
	DATA();
	virtual ~DATA();
	//! The method to do initialization
	/*!
	 *  Initialize particles and some settings
	 */
	void Initialization();
	//! The method to build neighbour list and coefficient list
	/*!
	 *  This function will call functions from octree search class and coefficient calculation class
	 */
	void Buildup_neigh_list_and_ceoff_list(const vector<double> &xp, const vector<double> &yp, const vector<double> &zp, const double &distance, const int num_of_par);
	//! The method to clear neighbour list and coefficient list
	void Clear_neigh_list_and_ceoff_list();
	//! The method to calculate coefficients for all particles
	/*!
	 *  \param neigh input The neighbour list of corresponding particle 
	 *  \param coeff1 input/output The coefficients for laplacian
	 *  \param coeff2 input/output the coefficients for dudx
	 *  \param coeff3 input/output the coefficients for dudy
	 *  \param coeff3 input/output the coefficients for dudz
	 */
	void GetLSCoefficient(const vector<int> &neigh, vector<double> &coeff1, vector<double> &coeff2, vector<double> &coeff3, vector<double> &coeff4);
	
public:
	//! Utilities
	//! The method to calculate cross product rslt = vec1 X vec2
	/*!
	 *  \param vec1 input The first vector
	 *  \param vec2 input The second vector
	 *  \param rslt output The result
	 */
	void Cross_Product(const vector<double> &vec1, const vector<double> &vec2, vector<double> &rslt);
	//! The method to get magnetic field
	/*!
	 *  /param x input x position of given point
	 *  /param y input y position of given point
	 *  /param z input z position of given point
	 *  /param B output The corresponding magnetic field
	 */
	void Get_MagneticFiled(const double x, const double y, const double z, vector<double> &B);

	//! The method to fetch number of particles
	int Get_num_of_par() {return m_num_of_par;}
	int Get_num_of_boundary_par() {return m_num_of_boundary_par;}
	//! The mothod to fetch phi
	vector<double>& Get_Phi() {return m_phi;}
	//! The method to fetch J
	vector<double>& Get_Jx() {return m_Jx;}
	vector<double>& Get_Jy() {return m_Jx;}
	vector<double>& Get_Jz() {return m_Jx;}
	//! The mothods to fetch position
	vector<double>& Get_x() {return m_xp;}
	double Get_x(int index) {return m_xp[index];}
	vector<double>& Get_y() {return m_yp;}
	double Get_y(int index) {return m_yp[index];}
	vector<double>& Get_z() {return m_zp;}
	double Get_z(int index) {return m_zp[index];}
	//! The mothods to fetch velocity
	vector<double>& Get_u() {return m_up;}
	double Get_u(int index) {return m_up[index];}
	vector<double>& Get_v() {return m_vp;}
	double Get_v(int index) {return m_vp[index];}
	vector<double>& Get_w() {return m_wp;}
	double Get_w(int index) {return m_wp[index];}

	//! The method to get boundary flag
	int Get_Boundary_Flag(int index) {return m_Boundary_Flag[index];}
	//! The method to fetch neighbour_list
	vector<vector<int> >& Get_neighbour_list() {return m_neighbour_list;}
	//! The methods to fetch coefficients list
	vector<vector<double> >& Get_coefficient_laplacian() {return m_coefficient_laplacian;}
	vector<vector<double> >& Get_coefficient_laplacian_boundary() {return m_coefficient_laplacian_boundary;}
	vector<vector<double> >& Get_coefficient_dudx() {return m_coefficient_dudx;}
	vector<vector<double> >& Get_coefficient_dudx_boundary() {return m_coefficient_dudx_boundary;}
	vector<vector<double> >& Get_coefficient_dudy() {return m_coefficient_dudy;}
	vector<vector<double> >& Get_coefficient_dudy_boundary() {return m_coefficient_dudx_boundary;}
	vector<vector<double> >& Get_coefficient_dudz() {return m_coefficient_dudz;}
	vector<vector<double> >& Get_coefficient_dudz_boundary() {return m_coefficient_dudz_boundary;}
	//! The methods to fetch normal
	vector<double>& Get_normal_x() {return m_normal_x;}
	vector<double>& Get_normal_y() {return m_normal_y;}
	vector<double>& Get_normal_z() {return m_normal_z;}
	//! The mothods to fetch compressible related quantities
	vector<double>& Get_energy() {return m_energy;}
	vector<double>& Get_rho() {return m_rho;}
	vector<double>& Get_pressure() {return m_pressure;}



private:
	//! Particle number
	int m_num_of_par, m_num_of_boundary_par;
	//! initial shortest particle distance
	double m_distance;
	//! Geometry quantity:normal of all particle
	vector<double> m_normal_x, m_normal_y, m_normal_z;
	//! Geometry quantity: Boundary particle or not. If the particle is boundary particle Boundary_flag[i] == 1;
	vector<int> m_Boundary_Flag;
    //! Particle position and velocity, density, volume, mass
    vector<double> m_xp, m_yp, m_zp, m_up, m_vp, m_wp, m_rho, m_vol, m_mass;
    //! Compressible code related quantities, energy density, pressure, temperature
    vector<double> m_energy, m_pressure, m_temperature;
    //! MHD related quantities
    vector<double> m_Jx, m_Jy, m_Jz, m_phi;
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
    //! The coefficient list for boundary particle
    /*!
     *  Store the coefficients of dudx, dudy, dudz, and laplacian of all particles, each coefficient is for the particle in the corresponding place of neighbour_list
     */
    vector<vector<double> > m_coefficient_laplacian_boundary, m_coefficient_dudx_boundary, m_coefficient_dudy_boundary, m_coefficient_dudz_boundary;
    //! The pointer to octree class
    /*!
     *  For building neighbour list
     */
    Octree *m_octree;
};

} /* namespace std */
#endif /* DATA_H_ */
