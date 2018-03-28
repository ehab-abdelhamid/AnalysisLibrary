/*
 * Utils.cpp
 *
 *  Created on: 3 Oct 2017
 *      Author: eabdelha
 */

#include<cstdlib>
#include<sstream>
#include <algorithm>
#include "Utils.h"

std::vector<Object*> sample(std::vector<Object*> list, float sampleR) {

	std::vector<Object*> rList;

	for(std::vector<Object*>::iterator it = list.begin(); it != list.end(); ++it) {
		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		if(r<=sampleR)
			rList.push_back(*it);
	}

	return rList;
}

std::vector<Object*> sample(std::vector<Object*> list, int numSamples) {

	if(numSamples>=list.size())
	{
		return std::vector<Object*>(list);
	}

	float ratio = ((float)numSamples)/list.size();
	return sample(list, ratio);
}

double dot(Object* obj1, Object* obj2) {
	std::vector<double>::iterator iter1 = obj1->getBeginning();
	std::vector<double>::iterator iter2 = obj2->getBeginning();
	double value = 0;

	for(; iter1 != obj1->getEnd(); ++iter1, ++iter2) {
		value+=((*iter1)*(*iter2));
	}

	return value;
}

Object* subtract(Object* obj1, Object* obj2) {
	Object* obj = new Object();

	std::vector<double>::iterator iter1 = obj1->getBeginning();
	std::vector<double>::iterator iter2 = obj2->getBeginning();

	for(; iter1 != obj1->getEnd(); ++iter1, ++iter2) {
		obj->addValue((*iter1)-(*iter2));
	}

	return obj;
}

string intToString(int a)
{
	stringstream sstmGF;
	sstmGF << a;
	return sstmGF.str();
}

string doubleToString(double a)
{
	stringstream sstmGF;
	sstmGF << a;
	return sstmGF.str();
}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

int getPos(char* str, char c, int start)
{
	for(int i=start;i<strlen(str);i++)
	{
		if(str[i]==c)
			return i;
	}
	return -1;
}

bool operator<(const orderedLists& a, const orderedLists& b) {
	return (a.getCurrentValue() > b.getCurrentValue());
}
