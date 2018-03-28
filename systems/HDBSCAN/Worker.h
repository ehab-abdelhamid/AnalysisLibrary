/*
 * Worker.h
 *
 *  Created on: 7 Feb 2018
 *      Author: eabdelha
 */

#ifndef WORKER_H_
#define WORKER_H_

#include "RandomBlocksSystem.h"

#define BUFFERSIZE 2000

class Worker {
private:
	RandomBlocksSystem rbs;
	static std::pair<int, int>* parseTaskMSG(char* str_message);
	static int nCorePoints;

public:
	int machineID;
	int nThreads;
	void start(int rank);
	static void *processing(Worker* worker, int id, RandomBlocksSystem* rbs);
};


#endif /* WORKER_H_ */
