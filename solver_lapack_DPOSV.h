/** 
 * \file	solver_lapack_ls.h
 * \brief	A wrap of dgelsd.f from lapack for solving the least square solution
 *  of Ax=b.
 *
 */
#ifndef  _LAPACK_DPOSV_H
#define  _LAPACK_DPOSV_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 *
 * \brief	A wrap of dgesvx.f from lapack for solving Ax=b.
 */
class LAPACK_DPOSV{
public:
	LAPACK_DPOSV();
	LAPACK_DPOSV(int n);
	~LAPACK_DPOSV();
	void Create(int n);

	void Set_A(int i, int j, double val);		// A[i][j]=val;
	void Add_A(int i, int j, double val);		// A[i][j]=A[i][j]+val;
	//void Set_x(int i, double val);			// x[i]=val;
	//void Set_x(double *p);					// x[i]=p[i];
	//void Add_x(int i, double val);			// x[i]=x[i]+val;
	void Get_x(double *p);						// get the x from ij_x to p.
	//void Get_x(double *p, int n, int *global_index);
	void Set_b(int i, double val);				// b[i]=val;
	void Set_b(double *b);
	void Add_b(int i, double val);				// b[i]=b[i]+val;

	//void SetMaxIter(int val){};
//	void GetFinalRelativeResidualNorm(double *rel_resid_norm);		// this error norm is not a relative residual norm
	//void GetNumIterations(int *num_iterations){};

	bool Solve(void);
	//virtual void Solve_withPureNeumann(void){};
	//virtual void Read_A(char *filename){};
	void Print_A(const char *filename);
	//virtual void Read_b(char *filename){};

	void Print_b(const char *filename);
	//virtual void Read_x(char *filename){};
	void Print_x(const char *filename);
	//virtual void test(void){};

private:
	// scalar arguments
	char m_UPLO;	        // input
	int m_N;		// input
	int m_NRHS;		// input
	double *m_A;	        // input/output
	int m_LDA;		// input
	double *m_B;	        // input/output
	int m_LDB;		// input
	int m_INFO;		// output

};

#endif  // #ifndef CLASS_SOLVER_H

