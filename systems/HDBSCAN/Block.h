/*
 * Block.h
 *
 *  Created on: 3 Oct 2017
 *      Author: eabdelha
 */

#ifndef BLOCK_H_
#define BLOCK_H_

#include<vector>
#include<map>
#include<set>
#include "Object.h"

class BlocksSystem;

using namespace std;

/**
 * a block of objects resulted from a tre-based partitioning technique
 */
class Block {
private:
	int marginLength;//the distance add for the boundary with neighbours
	int difaultNPoints = 3;
	map<int, Object*> elements;
	map<int, Object*> borderElems;
	std::pair<Object*, Object*>* projectionLine = NULL;
	Object* centerPoint = NULL;	//point in the middle for all points in the block
	double radius = -1;	//the distance from the centre point and the farthest point within this block
	void createProjectionLine();	//Create a projection line using elements belonging to this block
	Object* projectOnProjectionline(Object* );
	Block* left;
	Block* right;

public:
	Block();
	void addElemets(vector<Object*>, bool coreElement);
	void addElement(Object* element, bool coreElement);
	Object* getElementByID(int );
	Object* getProjection(Object* , std::pair<Object*, Object*>* );
	Object* pickRandomElement();
	Object* getFarthestPoint(Object* obj);
	int classify(Object*);
	void partition(vector<Block*>* smallestBlocks = NULL);
	int getNumberOfElements() { return elements.size(); }
	void calcCentreAndRadius();
	void computeCoreDistances(BlocksSystem*);
	void getCoreDistance(Object* , int , vector<pair<Object*,double>* >* );
	void print(ostream& out, bool printElems = false);
	int getMarginLength() { return marginLength; }
	set<Object* >* getObjectsSet();
	double getRadius() { return radius; }
	Object* getCenterPoint() { return centerPoint; }
	~Block();

};


#endif /* BLOCK_H_ */
