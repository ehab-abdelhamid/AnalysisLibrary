/*
 * RandomBlocksSystem.h
 *
 *  Created on: 9 Oct 2017
 *      Author: eabdelha
 */

#ifndef RANDOMBLOCKSSYSTEM_H_
#define RANDOMBLOCKSSYSTEM_H_

#include<vector>
#include<list>
#include "RandomBlock.h"
#include "CondenseHierarchy.h"
#include "GraphX.h"
#include "Clusters.h"

using namespace std;

class RandomBlocksSystem {
private:
	vector<RandomBlock* > blocks;
	vector<Object*> allObjects;
	void buildMST1(GraphX* , list<PossibleEdge*>* );
	void buildMSTApprox1(GraphX* , list<PossibleEdge*>* , double);
	list<PossibleEdge*>* buildLocalMSTofPairs(RandomBlock* b1, RandomBlock* b2);
	list<PossibleEdge*>* buildBipartiteMSTofPairs(RandomBlock* b1, RandomBlock* b2);

public:
	vector<Object*> allElements;
	void addBlock(RandomBlock* );
	RandomBlock* getBlock(int index);
	Object* getObject(int index);
	void computeCoreDistances(int );
	vector<RandomBlock*>::iterator getBegin() { return blocks.begin(); }
	vector<RandomBlock*>::iterator getEnd() { return blocks.end(); }
	RandomBlock* getBlockOfObject(Object* obj);
	void createBlocks(int );
	void assignElements(vector<Object*> allElements);
	void prepareBlocks();
	void prepareAllElements();
	void prepareSortedElements4Blocks();
	void buildMST();
	void buildMSTFromListOfEdges(vector<Object*> objectsList, GraphX* mst, list<PossibleEdge*>* total, bool needMap);
	void buildMSTFromListofPairs(GraphX* mst, list<PossibleEdge*>* total, vector<std::pair<int, int>*> blocksPairs);
	vector<PossibleEdge*>* buildBipartiteMSTofPairs_Fast(RandomBlock* b1, RandomBlock* b2);
	//void generateClusters(CH_Node* ch, double maxArea);
	Clusters* cluster(int minClusterSize);
	Clusters* clusterApprox1(int minClusterSize, double );
	Clusters* doClustering(GraphX* mst, list<PossibleEdge*>* edgesSorted, int minClusterSize);
	void print(ostream& out);
	~RandomBlocksSystem();
};


#endif /* RANDOMBLOCKSSYSTEM_H_ */
