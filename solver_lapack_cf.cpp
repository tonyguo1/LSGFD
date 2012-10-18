/**
 * \file	solver_lapack_ls.c
 * \brief	see also solver_lapack.h.
 *
 * b should always be a vector.
 */
#include "solver_lapack_cf.h"
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

  void FORTRAN_NAME(dpotrf)(char*,int*,double*,int*,int*);
#if defined(__cplusplus)
}
#endif


LAPACK_CF::LAPACK_CF()
{
	m_A = NULL;
}
LAPACK_CF::~LAPACK_CF()
{
	if(m_A!=NULL)
	{
		delete [] m_A;
	}
}
LAPACK_CF::LAPACK_CF(int m)
{
	Create(m);
}

void LAPACK_CF::Create(int m)
{
    m_UPLO = 'L';
	m_N = m;
	m_LDA = m;
	m_A = new double[m_LDA*m_N];

	int i;
	for(i=0; i<m_N*m_N; i++)
		m_A[i] = 0;
}

void LAPACK_CF::Set_A(int i, int j, double val)
{
	m_A[i+j*m_N] = val;
}
void LAPACK_CF::Add_A(int i, int j, double val)
{
	m_A[i+j*m_N] += val;
}
double LAPACK_CF::Get_L(int i,int j)
{
        return m_A[i+j*m_N];
}
//// this error norm is not a relative residual nor
//void LAPACK_CF::GetFinalRelativeResidualNorm(double *rel_resid_norm)
//{
//}

int LAPACK_CF::Solve(void)
{
        FORTRAN_NAME(dpotrf)(&m_UPLO,&m_N,m_A,&m_LDA,&m_INFO);

	//assert(m_INFO==0);
	return m_INFO;
}
