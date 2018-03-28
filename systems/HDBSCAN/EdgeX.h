/*
 * EdgeX.h
 *
 *  Created on: Mar 16, 2013
 *      Author: ehab
 *  This is for edge interface
 */

#ifndef EDGEX_H_
#define EDGEX_H_

#include "NodeX.h"

class EdgeX
{
private:
	double label;
	NodeX* otherNode;

public:
	EdgeX(double label, NodeX* otherNode);
	double getLabel() {return label;}
	NodeX* getOtherNode(){return otherNode;}
};

#endif /* EDGEX_H_ */
