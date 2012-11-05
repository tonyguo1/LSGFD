/*
 * EOSBASE.h
 *
 *  Created on: Nov 5, 2012
 *      Author: tongfei
 */

#ifndef EOSBASE_H_
#define EOSBASE_H_
class EOS_BASE {
public:
	EOS_BASE();
	virtual ~EOS_BASE();
	//method
	virtual void eos(const double&, const double&, double&, double&)=0;
	virtual void ThermalEnergy(const double&, const double&, double&)=0;
};

#endif /* EOSBASE_H_ */
