/**
 * \file	solver_lapack_ls.h
 * \brief	A wrap of dpotrf.f from lapack for Cholesky Factorization
 *
 */
#ifndef  _LAPACK_CF_H
#define  _LAPACK_CF_H

#include <stdio.h>
#include <stdlib.h>

/**
 *
 * \brief	A wrap of dpotrf.f from lapack for Cholesky Factorization.
 */
class LAPACK_CF {
public:
	LAPACK_CF();
	LAPACK_CF(int m);
	~LAPACK_CF();
	void Create(int ilower);

	void Set_A(int i, int j, double val);		// A[i][j]=val;
	void Add_A(int i, int j, double val);		// A[i][j]=A[i][j]+val;
	double Get_L(int i, int j);		// A[i][j]=A[i][j]+val;


	int Solve(void);

private:
	// scalar arguments
	char m_UPLO;	// input
	int m_N;		// input
	double *m_A;	// input/output
	int m_LDA;		// input
	int m_INFO;		// input

};

#endif  // #ifndef CLASS_SOLVER_H

