/*
 * EOSSPOLY.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "EOSSPOLY.h"
#include "math.h"

EOS_SPOLY::EOS_SPOLY(){
	Set_gamma(5.7);
	Set_pinf(50000.0);
	Set_einf(0.0);
}

EOS_SPOLY::~EOS_SPOLY() {
	// TODO Auto-generated destructor stub
}

void EOS_SPOLY::eos(const double &rho, const double &TE,double &p,double &cs)
{
	p = (m_gamma - 1.0)*(TE + m_einf) * rho - m_gamma * m_pinf;
	cs = sqrt(m_gamma*(p + m_pinf)/rho);
	return;
}

void EOS_SPOLY::ThermalEnergy(const double &rho, const double &p, double &TE)
{
	TE = (p + m_gamma * m_pinf) / (m_gamma - 1.0)/rho - m_einf;
	return;
}
