/*
 * Object.h
 *
 *  Created on: 27 Sep 2017
 *      Author: eabdelha
 */

#ifndef OBJECT_H_
#define OBJECT_H_

#include <iostream>
#include<vector>
#include "ANN/ANN.h"

class BlocksSystem;

using namespace std;

class Block;
class RandomBlock;

class Object {
public:
	Object() {id = statID; statID++;}
	double const& operator[](std::size_t index) const
	{
		return values[index];
	}
	void addValue(double a);
	void print(ostream& out);
	void printValues(ostream& out);
	double distanceTo(Object* , double maxDistance=-1);
	double coreDistanceTo(Object* , double maxDist=-1);
	std::vector<double>::iterator getBeginning() { return values.begin(); }
	std::vector<double>::iterator getEnd() { return values.end(); }
	double getValueByIndex(int i);
	double computelength();
	void addValuesFrom(Object* );
	void setValuesFrom(Object* );
	void zeroValues(int n);
	void divideValuesBy(float v);
	void calculateCoreDistance(Block* , BlocksSystem* bs, int);
	void calculateCoreDistance(RandomBlock* block, int nPoints);
	double getCoreDistance() { return this->coreDist; }
	void setCoreDistance(double value) { this->coreDist = value; }
	pair<Object*, double>* getKNNObject(Block* , BlocksSystem* bs, int);
	pair<Object*, double>* getKShortestRechability(Block* block, BlocksSystem* bs, int coreK, int k);
	vector<pair<Object*, double>*>* getObjectsWithDistance(Block* block, BlocksSystem* bs, double min, double max);
	long getID() { return id; }
	void setID(long id) { this->id = id; }
	static long getMaxID() { return statID; }
	void setBlock(Block* b) { block = b; }
	Block* getBlock() { return block; }
	int getValuesSize() { return values.size(); }
	ANNpoint getValuesAsANNPoint() { return &values[0]; }

private:
	long id;
	static long statID;
	vector<double> values;	//we assume the values for all objects follow the same ordering here
	double coreDist = -1;	//core distance of this point, -1 means it is not calculated
	Block* block = NULL;	//the block that this object belongs to as a core
};

#endif /* OBJECT_H_ */
