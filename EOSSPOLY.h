/*
 * EOSSPOLY.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef EOSSPOLY_H_
#define EOSSPOLY_H_

#include "EOSBASE.h"

class EOS_SPOLY: public EOS_BASE {
public:
	EOS_SPOLY();
	virtual ~EOS_SPOLY();
	virtual void eos(const double&, const double&, double&, double&);
	virtual void ThermalEnergy(const double&, const double&, double&);
	void Set_gamma(double gamma){m_gamma = gamma;}
	void Set_pinf(double pinf){m_pinf = pinf;}
	void Set_einf(double einf){m_einf = einf;}
	//data
    double m_gamma;
    double m_pinf;
    double m_einf;
};

#endif /* EOSSPOLY_H_ */
