/*
 * Master.h
 *
 *  Created on: 7 Feb 2018
 *      Author: eabdelha
 */

#ifndef MASTER_H_
#define MASTER_H_

#include <string>
#include "RandomBlocksSystem.h"
#include "PossibleEdge.h"

using namespace std;

class Master {
private:
	int nThreads;
	RandomBlocksSystem rbs;
	list<PossibleEdge*>* parseEdgesMSG(char* str_message);
	void sendTask(std::pair<int, int>* task, int currentWorker);

public:
        string output;
	void start(string fileName, bool ignoreFirstLine, int numMachine, int nThreads, bool approximate, int numPartitions, int nCorePoints, int minClusterSize);
};



#endif /* MASTER_H_ */
