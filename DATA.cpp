/*
 * DATA.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: tongfei
 */

#include "DATA.h"
#include "octree.h"
#include <math.h>
namespace std {

DATA::DATA() {
	m_octree = NULL;
	m_num_of_par = 0;
	m_num_of_boundary_par = 0;
	m_distance = 0;
	m_max_step = 0;
	m_max_time = 0;
	// TODO Auto-generated constructor stub

}

DATA::~DATA() {
	// TODO Auto-generated destructor stub
}

void DATA::Initialization(){
	/** For testing run, hard coded initialization**/
	//! Lower bound and upper bound
	double L[3] = {-0.5, -0.5, -8};
	double U[3] = {0.5, 0.5, 0};
	double cen[3] = {0, 0, -4}, R = 0.5, length = 3.5, dist = 0;
	m_distance = 0.1; 
	double xr, yr, zr;
	int num_par_x, num_par_y, num_par_z;
	num_par_x = static_cast<int>((U[0] - L[0])/m_distance + 1e-7) + 1;
	num_par_y = static_cast<int>((U[1] - L[1])/m_distance + 1e-7) + 1;
	num_par_z = static_cast<int>((U[2] - L[2])/m_distance + 1e-7) + 1;
	for (int i = 0; i < num_par_x; i++)
		for (int j = 0; j < num_par_y; j++)
			for (int k = 0; k < num_par_z; k++){
				double x = L[0] + i * m_distance, y = L[1] + j * m_distance, z = L[2] + k * m_distance;
				xr = x - cen[0];
				yr = y - cen[1];
				zr = z - cen[2];
				if (zr > length){
					dist = xr * xr + yr * yr + (zr - length) * (zr - length);
					dist = sqrt(dist) - R;
				}
				else if (zr > -length){
					dist = xr * xr + yr * yr;
					dist = sqrt(dist) - R;
				}
				else{
					dist = xr * xr + yr * yr + (zr + length) * (zr + length);
					dist = sqrt(dist) - R;
				}
				if (dist <= 0){
					//! Inside computational domain
					m_xp.push_back(x);
					m_yp.push_back(y);
					m_zp.push_back(z);
					m_up.push_back(0.21);
					m_vp.push_back(0.21);
					m_wp.push_back(0.21);
					m_rho.push_back(13);
					m_pressure.push_back(1000);
					//m_energy.push_back(eos->energy(13,1000));
					m_num_of_par++;
				}
			}
}

void DATA::Buildup_neigh_list_and_ceoff_list(const vector<double> &xp, const vector<double> &yp, const vector<double> &zp, const double &distance, const int num_of_par){
	int N = num_of_par;
	/** 5 is optimum for search */
	m_octree = new Octree(xp, yp, zp, 5, N);
	vector<double> coeff;
	vector<double> coeff1(30,0);
	vector<double> coeff2(30,0);
	vector<double> coeff3(30,0);
	vector<double> coeff4(30,0);
	vector<int> neigh;
	m_neighbour_list.assign(num_of_par,neigh);
	m_coefficient_laplacian.assign(num_of_par,coeff);
	m_coefficient_dudx.assign(num_of_par,coeff);
	m_coefficient_dudy.assign(num_of_par,coeff);
	m_coefficient_dudz.assign(num_of_par,coeff);
	m_Boundary_Flag.assign(num_of_par,0);
	m_normal_x.assign(num_of_par,0);
	m_normal_y.assign(num_of_par,0);
	m_normal_z.assign(num_of_par,0);
	for (int i_index = 0; i_index < num_of_par; i_index++){
		/** Get the neighbour list*/
		vector<int> octree_search_result;
		vector<double> octree_search_distance;
		int search_num = m_octree->searchNeighbor(xp[i_index], yp[i_index], zp[i_index], distance, &octree_search_result, &octree_search_distance );
		pair<double, int> pa;
		vector<pair<double, int> > vecs;
		vector<int>::iterator it_r = octree_search_result.begin();
		vector<double>::iterator it_d = octree_search_distance.end();
		for (int i = 0; i < search_num; i++){
			pa.first = *it_d;
			pa.second = *it_r;
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


		/** Get the coefficients*/
		GetLSCoefficient(neigh, coeff1, coeff2,coeff3,coeff4);
		neigh.clear();
	}
}

void DATA::Clear_neigh_list_and_ceoff_list(){
	int size = m_coefficient_laplacian.size();
	if (size != 0){
		for (int i = 0; i < size; i++){
			m_coefficient_laplacian[i].clear();
			m_coefficient_dudx[i].clear();
			m_coefficient_dudy[i].clear();
			m_coefficient_dudz[i].clear();
		}
		m_coefficient_laplacian_boundary.clear();
		m_coefficient_dudx_boundary.clear();
		m_coefficient_dudy_boundary.clear();
		m_coefficient_dudz_boundary.clear();
	}
	size = m_coefficient_laplacian_boundary.size();
	if (size != 0){
			for (int i = 0; i < size; i++){
				m_coefficient_laplacian_boundary[i].clear();
				m_coefficient_dudx_boundary[i].clear();
				m_coefficient_dudy_boundary[i].clear();
				m_coefficient_dudz_boundary[i].clear();
			}
			m_coefficient_laplacian_boundary.clear();
			m_coefficient_dudx_boundary.clear();
			m_coefficient_dudy_boundary.clear();
			m_coefficient_dudz_boundary.clear();
		}
}

void DATA::GetLSCoefficient(const vector<int> &neigh, vector<double> &coeff1, vector<double> &coeff2, vector<double> &coeff3, vector<double> &coeff4){
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
		g[i] = m_zp[neigh[i]] - m_zp[neigh[0]];
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
	double angle = *(min_element(angles+1, angles+n));
	m_angle.push_back(angle);
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
	m_coefficient_laplacian[neigh[0]].assign(coeff1.begin(),coeff1.end());
	m_coefficient_dudx[neigh[0]].assign(coeff2.begin(),coeff2.end());
	m_coefficient_dudy[neigh[0]].assign(coeff3.begin(),coeff3.end());
	m_coefficient_dudz[neigh[0]].assign(coeff4.begin(),coeff4.end());

	// Calculating coefficients with ghost particle
	if (angle > PI / 4){
		m_Boundary_Flag[neigh[0]] = 1;
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

const char* DATA::right_flush(
		int	   n,
		int	   ndigits)
{
	static	char	s[20];
	int		i;

	if (n == 0)
		ndigits--;
	for (i = n; i > 0; i /= 10) ndigits--;

	s[0] = '\0';
	for (;ndigits > 0; ndigits--)
		(void) sprintf(s,"%s0",s);
	(void) sprintf(s,"%s%d",s,n);
	return s;
}

void DATA::Print(const double t, const int step, const char* outputname){
	size_t i;
	char filename[200];
	FILE *outfile;
	sprintf(filename,"vtkoutput_%s",right_flush(step,7));
	sprintf(filename,"%s_%s",filename,outputname);
	//if (PN * PM * PL != 1)
	//	sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.vtk",filename);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.8f\n",t);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	fprintf(outfile,"POINTS %lu double\n",m_num_of_par);
	for (i = 0; i < m_num_of_par; i++)
		fprintf(outfile,"%.16g %.16g %.16g\n",m_xp[i],m_yp[i],m_zp[i]);
	fprintf(outfile,"POINT_DATA %d\n",m_num_of_par);
	fprintf(outfile,"VECTORS Current double\n");
	for (i = 0; i < m_num_of_par; i++)
		fprintf(outfile,"%.16g %.16g %.16g\n",m_up[i],m_vp[i],m_wp[i]);
	fclose(outfile);
}

void DATA::Cross_Product(const vector<double> &vec1, const vector<double> &vec2, vector<double> &rslt){
	rslt[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	rslt[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	rslt[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void DATA::Get_MagneticFiled(const double x, const double y, const double z, vector<double> &B){
	B[0] = 14.1 * sqrt(0.5*(1.0 - tanh((z - 1.5)/0.62)));
	B[1] = 0;
	B[2] = 0;
}

double DATA::Get_dt(){
	return 0.005;
}

void DATA::Print_angle(){
	FILE *fp;
	fp = fopen("angle.dat","w+");
	for (int i = 0; i < m_num_of_par; i++)
		fprintf(fp,"%d %lf\n", i, m_angle[i]);
	fclose(fp);
}
} /* namespace std */
