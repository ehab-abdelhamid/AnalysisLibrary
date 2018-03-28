/*
 * RandomBlock.h
 *
 *  Created on: 27 Nov 2017
 *      Author: eabdelha
 */

#ifndef RANDOMBLOCK_H_
#define RANDOMBLOCK_H_

#include<vector>
#include<list>
#include <boost/pending/disjoint_sets.hpp>
#include "Object.h"
#include "ClosePartition.h"

using namespace std;

class RandomBlock {
private:
	int id;
	vector<Object*> elements;//elements assigned to this partition
	vector<Object*> sortedElements;//elements assigned to this partition but sorted based on core distance
	vector<Object*> objects;//all available data from the system
//	vector<ClosePartition*> partitions;
//	int defaultNPoints = 3;
	ANNkd_tree*	kdTreeObjects = NULL;
	ANNpointArray dataPts1;
	ANNkd_tree*	kdTreeElements = NULL;
	ANNpointArray dataPts2;

public:
	RandomBlock();
	~RandomBlock();
	void addElement(Object* obj);
	void addObject(Object* obj);
//	void partition(long r);
	void prepare();
//	vector<PossibleEdge*> getPairsWithDistBound(double min, double max);
	void print(ostream& out);
	void computeCoreDistances(int );
	void prepareSortedElements();
	void getCoreDistance(Object* , int , vector<pair<Object*,double>* >* );
	vector<Object*> getElements() { return elements; }
	void buildLocalMST();
	vector<PossibleEdge*>* buildLocalMST_Prims();
	vector<Object*>::iterator getElementsIter() { return elements.begin(); }
	vector<Object*>::iterator getElementsEndIter() { return elements.end(); }
	int getNumberofElements() { return elements.size(); }
	std::pair<Object*, double>* getClosestElement_CoreDistance(Object* obj);
	void setID(int id) { this->id = id; }
	int getID() { return id; }

	static long time1;
};


#endif /* RANDOMBLOCK_H_ */
