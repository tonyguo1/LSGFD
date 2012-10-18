/*
 * DATA.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "DATA.h"
#include "octree.h"

namespace std {

DATA::DATA() {
	// TODO Auto-generated constructor stub

}

DATA::~DATA() {
	// TODO Auto-generated destructor stub
}

void DATA::Buildup_neigh_list_and_ceoff_list(){
	int N = m_num_of_par;
	/** 5 is optimum for search */
	m_octree = new Octree(m_xp, m_yp, m_zp, N, 5);
	vector<double> coeff1(30,0);
	vector<double> coeff2(30,0);
	vector<double> coeff3(30,0);
	vector<double> coeff4(30,0);
	vector<int> neigh;
	m_neighbour_list.assign(m_num_of_par,neigh);
	m_coefficient_laplacian.assign(m_num_of_par,coeff1);
	m_coefficient_dudx.assign(m_num_of_par,coeff1);
	m_coefficient_dudx.assign(m_num_of_par,coeff1);
	m_coefficient_dudx.assign(m_num_of_par,coeff1);
	Boundary_Flag.assign(m_num_of_par,0);
	m_normal_x.assign(m_num_of_par,0);
	m_normal_y.assign(m_num_of_par,0);
	m_normal_z.assign(m_num_of_par,0);
	for (int i_index = 0; i_index < m_num_of_par; i_index++){
		/** Get the neighbour list*/
		vector<int> *octree_search_result = new vector<int>;
		vector<double> *octree_search_distance = new vector<int>;
		int search_num = m_octree->searchNeighbor(m_xp[i_index], m_yp[i_index], m_zp[i_index], 2 * m_distance, octree_search_result, octree_search_distance );
		neigh.push_back(i_index);
		pair<double, int> pa;
		vector<pair<double, int> > vecs;
		vector<int>::iterator it_r = octree_search_result->begin();
		vector<double>::iterator it_d = octree_search_distance->end();
		for (int i = 0; i < search_num; i++){
			pa.first = *it_r;
			pa.second = *it_d;
			vecs.push_back(pa);
			it_r++;
			it_d++;
		}
		vector<pair<double, int> >::iterator it_v = vecs.begin();;
		sort(vecs.begin(),vecs.end());
		int num = 1;
		while (num < 28 && it_v != vecs.end()){
			neigh.push_back((*it_v).second);
			num++;
			++it_v;
		}
		m_neighbour_list[i_index] = neigh;
		neigh.clear();
		delete octree_search_result;
		delete octree_search_distance;


		/** Get the coefficients*/
		GetLSCoefficient(neigh, coeff1, coeff2,coeff3,coeff4);
	}
}

void DATA::clear_neigh_list_and_ceoff_list(){
	int size = m_coefficient_laplacian_boundary.size();
	if (size != 0){
			for (int i = 0; i < size; i++){
				m_coefficient_laplacian_boundary[i].clear();
				m_coefficient_dudx_boundary[i].clear();
				m_coefficient_dudx_boundary[i].clear();
				m_coefficient_dudx_boundary[i].clear();
			}
			m_coefficient_laplacian_boundary.clear();
			m_coefficient_dudx_boundary.clear();
			m_coefficient_dudx_boundary.clear();
			m_coefficient_dudx_boundary.clear();
		}
}

void DATA::GetLSCoefficient(vector<int> neigh, vector<double> coeff1, vector<double> coeff2,vector<double> coeff3, vector<double> coeff4){
	int n = neigh.size() - 1;
	double normal[3] = {0};
	double h[30] = {0};
	double k[30] = {0};
	double g[30] = {0};
	for ( int i = 1; i <= n; i++ ){
		h[i] = m_xp[neigh[i]] - m_xp[neigh[0]];
		normal[0] -= h[i];
		k[i] = m_yp[neigh[i]] - m_yp[neigh[0]];
		normal[1] -= k[i];
		k[i] = m_zp[neigh[i]] - m_zp[neigh[0]];
		normal[2] -= g[i];
	}
	m_normal_x[neigh[0]] = normal[0];
	m_normal_y[neigh[0]] = normal[1];
	m_normal_z[neigh[0]] = normal[2];
	double angles[30] = {0};
	double length2 = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	for (int i = 1; i <= n; i++){
		double dot = (h[i] * normal[0] + k[i] * normal[1] + g[i] * normal[2]);
		double length1 = sqrt(h[i] * h[i] + k[i] * k[i] + g[i] * g[i]);
		if (length2 < 1e-15)
			angles[i] = 0;
		else{
			double costheta  = dot / (length1 * length2);
			angles[i] = acos(costheta);
		}
	}
	double angle = min_element(angles,angles+30);
	//initialize A
	double A[10][10];
	for (int i = 1; i <= 9; i++)
		for (int j = 1; j <= 9; j++)
			A[i][j] = 0;
	for (int index = 1; index <= n; index++){
		double a[10] = {0,h[index],k[index],g[index],0.5*h[index]*h[index],0.5*k[index]*k[index],0.5*g[index]*g[index],
				h[index]*k[index],h[index]*g[index],k[index]*g[index]};
		for (int i = 1; i <= 9; i++)
			for (int j = 1; j <= 9; j++)
				A[i][j] += a[i] * a[j];
	}



	//get L and L^T
	LAPACK_CF lapack_cf(9);
	for (int i = 1; i <= 9; i++)
		for (int j = 1; j <= 9; j++)
			lapack_cf.Set_A(i-1,j-1,A[i][j]);
	int status = lapack_cf.Solve();

	//get value of L
	double l[10][10];
	for (int i = 1; i < 10; i++)
		for (int j = 1; j <=i; j++)
			l[i][j] = lapack_cf.Get_L(i-1,j-1);
	delete lapack_cf;
	//build M
	double M[10][10];
	double c[10] = {0};
	double d[30][10];
	for (int index = 1; index <= n; index++){
		double a[10] = {0,h[index],k[index],g[index],0.5*h[index]*h[index],0.5*k[index]*k[index],0.5*g[index]*g[index],
				h[index]*k[index],h[index]*g[index],k[index]*g[index]};
		for (int i = 1; i <=9; i++){
			d[index][i]=a[i];
			c[i]+=d[index][i];
		}
	}
	for ( int i = 1; i <= 9; i++ )
		for ( int j = 9; j >= 1; j-- ){
			if (j == i)
				M[i][j] = 1 / l[i][i];
			else if (j < i){
				double sum = 0;
				for (int k = j; k <= i - 1; k++){
					sum += l[i][k] * M[k][j];
				}
				M[i][j] = -1 * sum / l[i][i];
			}
			else
				M[i][j] = 0;
		}
	double MT[10][10], AI[10][10];
	for (int i = 1; i< 10; i++)
		for (int j = 1; j < 10; j++){
			MT[i][j] = M[j][i];
			AI[i][j] = 0;
		}
	for (int i = 1; i < 10; i++)
		for (int j = 1; j < 10; j++)
			for (int kk = 1; kk < 10; kk++)
				AI[i][j] += MT[i][kk] * M[kk][j];

	//calculate coefficients
	for (int i = 1; i <= 9; i++){
		coeff1[0] += -(AI[4][i] + AI[5][i] + AI[6][i])*c[i];
		coeff2[0] += -AI[1][i]*c[i];
		coeff3[0] += -AI[2][i]*c[i];
		coeff4[0] += -AI[3][i]*c[i];
	}
	for ( int j = 1; j <= n; j++){
		for ( int i = 1; i <= 9; i++ ){
			coeff1[j] += (AI[4][i] + AI[5][i] + AI[6][i])* d[j][i];
			coeff2[j] += AI[1][i] * d[j][i];
			coeff3[j] += AI[2][i] * d[j][i];
			coeff4[j] += AI[3][i] * d[j][i];
		}
	}
	m_coefficient_laplacian[neigh[0]] = coeff1;
	m_coefficient_dudx[neigh[0]] = coeff2;
	m_coefficient_dudy[neigh[0]] = coeff3;
	m_coefficient_dudz[neigh[0]] = coeff4;

	// Calculating coefficients with ghost particle
	if (angles < PI / 4){
		Boundary_Flag[neigh[0]] = 1;
		n++;
		h[n] = normal[0] / (n-1);
		k[n] = normal[1] / (n-1);
		g[n] = normal[2] / (n-1);

		//initialize A
		for (int i = 1; i <= 9; i++)
			for (int j = 1; j <= 9; j++)
				A[i][j] = 0;
		for (int index = 1; index <= n; index++){
			double a[10] = {0,h[index],k[index],g[index],0.5*h[index]*h[index],0.5*k[index]*k[index],0.5*g[index]*g[index],
					h[index]*k[index],h[index]*g[index],k[index]*g[index]};
			for (int i = 1; i <= 9; i++)
				for (int j = 1; j <= 9; j++)
					A[i][j] += a[i] * a[j];
		}

		//get L and L^T
		for (int i = 1; i <= 9; i++)
			for (int j = 1; j <= 9; j++)
				lapack_cf.Set_A(i-1,j-1,A[i][j]);
		int status = lapack_cf.Solve();

		//get value of L
		for (int i = 1; i < 10; i++)
			for (int j = 1; j <=i; j++)
				l[i][j] = lapack_cf.Get_L(i-1,j-1);
		delete lapack_cf;
		//build M
		for (int i = 1; i < 10; i++)
			c[i] = 0;
		for (int index = 1; index <= n; index++){
			double a[10] = {0,h[index],k[index],g[index],0.5*h[index]*h[index],0.5*k[index]*k[index],0.5*g[index]*g[index],
					h[index]*k[index],h[index]*g[index],k[index]*g[index]};
			for (int i = 1; i <=9; i++){
				d[index][i]=a[i];
				c[i]+=d[index][i];
			}
		}
		for ( int i = 1; i <= 9; i++ )
			for ( int j = 9; j >= 1; j-- ){
				if (j == i)
					M[i][j] = 1 / l[i][i];
				else if (j < i){
					double sum = 0;
					for (int k = j; k <= i - 1; k++){
						sum += l[i][k] * M[k][j];
					}
					M[i][j] = -1 * sum / l[i][i];
				}
				else
					M[i][j] = 0;
			}
		double MT[10][10], AI[10][10];
		for (int i = 1; i< 10; i++)
			for (int j = 1; j < 10; j++){
				MT[i][j] = M[j][i];
				AI[i][j] = 0;
			}
		for (int i = 1; i < 10; i++)
			for (int j = 1; j < 10; j++)
				for (int kk = 1; kk < 10; kk++)
					AI[i][j] += MT[i][kk] * M[kk][j];

		//calculate coefficients
		for (int i = 1; i <= 9; i++){
			coeff1[0] += -(AI[4][i] + AI[5][i] + AI[6][i])*c[i];
			coeff2[0] += -AI[1][i]*c[i];
			coeff3[0] += -AI[2][i]*c[i];
			coeff4[0] += -AI[3][i]*c[i];
		}
		for ( int j = 1; j <= n; j++){
			for ( int i = 1; i <= 9; i++ ){
				coeff1[j] += (AI[4][i] + AI[5][i] + AI[6][i])* d[j][i];
				coeff2[j] += AI[1][i] * d[j][i];
				coeff3[j] += AI[2][i] * d[j][i];
				coeff4[j] += AI[3][i] * d[j][i];
			}
		}
		m_coefficient_laplacian_boundary.push_back(coeff1);
		m_coefficient_dudx_boundary.push_back(coeff2);
		m_coefficient_dudy_boundary.push_back(coeff3);
		m_coefficient_dudz_boundary.push_back(coeff4);
		m_num_of_boundary_par++;
	}
}
} /* namespace std */
