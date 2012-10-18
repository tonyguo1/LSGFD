/**
 * \file	solver_lapack_DPOSV.c
 * \brief	see also solver_lapack.h.
 *
 * b should always be a vector.
 */
#include "solver_lapack_DPOSV.h"
#include <assert.h>
#include <algorithm>
// the following FROTRAN_NAME is copied from FronTier
#if defined(cray) || defined(_AIX)
#	define 	FORTRAN_NAME(a)		a
#else
#	define	FORTRAN_NAME(a)		a ## _
#endif

#if defined(__cplusplus)
extern "C" {
#endif

void FORTRAN_NAME(dposv)(char*,int*,int*,double*,int*,double*,int*,int*);

#if defined(__cplusplus)
}
#endif


LAPACK_DPOSV::LAPACK_DPOSV()
{
	m_A = NULL;
}
LAPACK_DPOSV::~LAPACK_DPOSV()
{
	if(m_A!=NULL)
	{
		delete [] m_A;
		delete [] m_B;
	}
}
LAPACK_DPOSV::LAPACK_DPOSV(int n)
{
	Create(n);
}

void LAPACK_DPOSV::Create(int n)
{
	m_UPLO = 'L';
	m_N = n;
	m_NRHS = 1;
	m_LDA = m_N;
	m_A = new double[m_LDA*m_N];
	m_LDB = m_N;
	m_B = new double[m_N];
	
	int i;
	for(i=0; i<m_LDB*m_NRHS; i++)
	{
		m_B[i] = 0;
	}
	for(i=0; i<m_N*m_N; i++)
		m_A[i] = 0;
}

void LAPACK_DPOSV::Set_A(int i, int j, double val)
{
	m_A[i+j*m_N] = val;
}
void LAPACK_DPOSV::Add_A(int i, int j, double val)
{
	m_A[i+j*m_N] += val;
}
void LAPACK_DPOSV::Get_x(double *x)
{
	int i;
	for(i=0; i<m_N; i++)
		x[i] = m_B[i];
}
void LAPACK_DPOSV::Set_b(int i, double val)
{
	m_B[i] = val;
}
void LAPACK_DPOSV::Set_b(double *b)
{
	int i;
	for(i=0; i<m_LDB; i++)
		m_B[i] = b[i];
}
void LAPACK_DPOSV::Add_b(int i, double val)
{
	m_B[i] += val;
}


bool LAPACK_DPOSV::Solve(void)
{
  //FORTRAN_NAME(dgels)(&m_TRANS,&m_M,&m_N,&m_NRHS,m_A,&m_LDA,m_B,&m_LDB,m_WORK,&m_LWORK,&m_INFO);
  FORTRAN_NAME(dposv)(&m_UPLO,&m_N,&m_NRHS,m_A,&m_LDA,m_B,&m_LDB,&m_INFO);

	assert(m_INFO==0);
	if(m_INFO!=0)
		return false;
	else
		return true;
}
void LAPACK_DPOSV::Print_A(const char *filename)
{
	int i, j;

	for(j=0; j<m_N; j++)
		for(i=0; i<m_N; i++)
			printf("A[%d][%d] = %f \n", i, j, m_A[i+j*m_N]);
}
void LAPACK_DPOSV::Print_b(const char *filename)
{
	int i;
	for(i=0; i<m_LDB; i++)
		printf("b[%d] = %f \n", i, m_B[i]);
}
void LAPACK_DPOSV::Print_x(const char *filename)
{
	int i;
	for(i=0; i<m_N; i++)
		printf("x[%d] = %f \n", i, m_B[i]);
}
