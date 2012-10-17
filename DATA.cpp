/*
 * DATA.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "DATA.h"

namespace std {

DATA::DATA() {
	// TODO Auto-generated constructor stub

}

DATA::~DATA() {
	// TODO Auto-generated destructor stub
}

void DATA::Buildup_neigh_list_and_ceoff(){
	int N = m_num_of_par;
	double *xp = new double[N];
	double *yp = new double[N];
	double *zp = new double[N];
	for (int i_index = 0; i_index < m_num_of_par; i_index++){
		xp[i_index] = m_xp[i_index];
		yp[i_index] = m_yp[i_index];
		zp[i_index] = m_zp[i_index];
	}
	/** 5 is optimum for search */
	DefaultCPUParticles particles(xp, yp, zp , N, 5);
	delete [] xp;
	delete [] yp;
	delete [] zp;
	for (size_t i_index = 0; i_index < m_num_of_par; i_index++){
		/** Get the neighbour list*/
		vector<int> *octree_search_result = new vector<int>;
		vector<double> *octree_search_distance = new vector<int>;
		int search_num = particles.searchNeighbor(m_xp[i_index], m_yp[i_index], m_zp[i_index], 2 * m_distance, (void **)(&octree_search_result), (void **)(&octree_search_distance) );
		size_t num = 1;
		vector<int> neigh;
		neigh.push_back(i_index);
		pair<double, int> pa;
		vector<pair<double, int> > vecs;
		vector<int>::iterator it_r = octree_search_result->begin();
		vector<double>::iterator it_d = octree_search_distance->end();
		for (int i = 0; i < search_num; i++){
			pair.first = *it_r;
			pair.second = *it_d;
			vecs.push_back(pa);
			it_r++;
			it_d++;
		}
		vector<pair<double, int> >::iterator it_v = vecs.begin();;
		sort(vecs.begin(),vecs.end());
		while (num < 28 && it_v != vecs.end()){
			neigh.push_back(*it_v.second);
			num++;
			++it_v;
		}
		neighbour_list.push_back(neigh);
		neigh.clear();
		delete octree_search_result2;
		delete octree_search_distance;


		/** Get the coefficients*/
		vector<double> coeff1(30,0);
		vector<double> coeff2(30,0);
		vector<double> coeff3(30,0);
		vector<double> coeff4(30,0);

	}
}
} /* namespace std */
