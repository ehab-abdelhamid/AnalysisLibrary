/*
 * Utils.h
 *
 *  Created on: 3 Oct 2017
 *      Author: eabdelha
 */

#ifndef UTILS_H_
#define UTILS_H_

#include<vector>
#include<list>
#include<queue>
#include "Object.h"
#include "PossibleEdge.h"

std::vector<Object*> sample(std::vector<Object*> list, float sampleR);
std::vector<Object*> sample(std::vector<Object*> list, int numSamples);
double dot(Object*, Object*);
Object* subtract(Object*, Object*);

string intToString(int a);
string doubleToString(double a);

int getPos(char* str, char c, int start = 0);

char* getCmdOption(char ** begin, char ** end, const std::string & option);

class orderedLists{
private:
	list<PossibleEdge*>* edges;
	list<PossibleEdge*>::iterator currentPosition;

public:
	orderedLists() {
		edges = NULL;
	}

	orderedLists(list<PossibleEdge*>* list) {
		edges = list;
		currentPosition = edges->begin();
	}

	PossibleEdge* getCurrentEdge() {
		return (*currentPosition);
	}

	bool advance() {
		currentPosition++;
		if(currentPosition!=edges->end())
			return true;
		else
			return false;
	}

	double getCurrentValue() const {
		return (*currentPosition)->dist;
	}

	int getSize() {
		return edges->size();
	}

	~orderedLists() {
	}
};

bool operator<(const orderedLists& a, const orderedLists& b);

#endif /* UTILS_H_ */
