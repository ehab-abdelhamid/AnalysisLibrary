/*
 * BlocksSystem.h
 *
 *  Created on: 9 Oct 2017
 *      Author: eabdelha
 */

#ifndef BLOCKSSYSTEM_H_
#define BLOCKSSYSTEM_H_

#include<vector>
#include "Block.h"
#include "CondenseHierarchy.h"

/////////////////////////////////////////
// I am not currently using this class //
/////////////////////////////////////////

using namespace std;

class BlocksSystem {
private:
	vector<Block* > blocks;

public:
	void addBlock(Block* );
	void computeCoreDistances();
	vector<Block*>::iterator getBegin() { return blocks.begin(); }
	vector<Block*>::iterator getEnd() { return blocks.end(); }
	Block* getBlockOfObject(Object* obj);
	void prepareBlocks();
	void buildMST();
	void generateClusters(CH_Node* ch, double maxArea);
	void print(ostream& out);
};


#endif /* BLOCKSSYSTEM_H_ */
