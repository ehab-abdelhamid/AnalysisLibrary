/*
 * PossibleEdge.cpp
 *
 *  Created on: 30 Nov 2017
 *      Author: eabdelha
 */

#include "PossibleEdge.h"
#include "Object.h"

PossibleEdge::PossibleEdge(Object* src, Object* dest, double dist) {
	this->src = src;
	this->dest = dest;
	this->dist = dist;
}

void PossibleEdge::print(ostream& out) {
	cout<<src->getID()<<", "<<dest->getID()<<", "<<this->dist<<endl;
}

void PossibleEdge::swapValues(PossibleEdge* a, PossibleEdge* b) {
	Object* tempSrc;
	Object* tempDest;
	double tempDist;
	tempSrc = a->src;tempDest = a->dest;tempDist = a->dist;
	a->src = b->src;a->dest = b->dest;a->dist = b->dist;
	b->src = tempSrc;b->dest = tempDest;b->dist = tempDist;
}

