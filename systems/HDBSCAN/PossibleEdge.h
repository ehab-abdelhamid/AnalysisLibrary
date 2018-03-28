/*
 * PossibleEdge.h
 *
 *  Created on: 30 Nov 2017
 *      Author: eabdelha
 */

#ifndef POSSIBLEEDGE_H_
#define POSSIBLEEDGE_H_

#include<iostream>

using namespace std;

class Object;

class PossibleEdge {
public:
	Object* src;	//source object
	Object* dest;	//destination object
	double dist;	//distance of the edge

	PossibleEdge(Object* src, Object* dest, double dist);
	void print(ostream& out);
	static void swapValues(PossibleEdge* a, PossibleEdge* b);
};


#endif /* POSSIBLEEDGE_H_ */
