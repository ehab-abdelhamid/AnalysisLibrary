/*
 * Master.cpp
 *
 *  Created on: 7 Feb 2018
 *      Author: eabdelha
 */

#include <mpi.h>
#include <vector>
#include <string.h>
#include <fstream>
#include "Master.h"
#include "Settings.h"
#include "CSVReader.h"
#include "Settings.h"
#include "RandomBlocksSystem.h"
#include "PossibleEdge.h"
#include "Utils.h"

using namespace std;

list<PossibleEdge*>* Master::parseEdgesMSG(char* str_message) {
	list<PossibleEdge*>* rList = new list<PossibleEdge*>();
	int pos = 1;
	char block1Str[10];
	char block2Str[10];
	char distStr[20];
	while(true) {
		//get the first number of the pair
		int pos1 = getPos(str_message, '|', pos+1);

		strncpy (block1Str, str_message+pos+1, pos1-pos);
		int obj1 = atoi(block1Str);

		pos = getPos(str_message, '|', pos1+1);
		strncpy (block2Str, str_message+pos1+1, pos-pos1);
		int obj2 = atoi(block2Str);

		pos1 = getPos(str_message, ',', pos+1);
		if(pos1==-1) {
//			if((strlen(str_message)-pos)>=20)
//				cout<<"size is larger1: "<<(strlen(str_message)-pos)<<endl<<flush;
			strncpy (distStr, str_message+pos+1, strlen(str_message)-pos);
		}
		else {
//			if((pos1-pos)>=20) {
//				str_message[pos+100]='\0';
//				cout<<"size is larger2: "<<(pos1-pos)<<" "<<(str_message+pos+1)<<endl<<flush;
//				cout<<"size is larger22: "<<(pos1-pos)<<" "<<(str_message+pos-100)<<endl<<flush;
//			}
			strncpy (distStr, str_message+pos+1, pos1-pos);
		}
		double dist = strtod(distStr, NULL);

		rList->push_back(new PossibleEdge(rbs.getObject(obj1), rbs.getObject(obj2), dist));

		pos = pos1;

		if(pos==-1)
			break;
	}

	return rList;
}

void Master::sendTask(std::pair<int, int>* task, int currentWorker) {
	std::stringstream msgStr;

	msgStr<<"t";
	msgStr<<","<<task->first<<"|"<<task->second<<"\n";

	const std::string& tmp = msgStr.str();
	const char* cstr = tmp.c_str();

	int cnt=strlen(cstr);//+1;
	MPI_Status status;
	int rank = (currentWorker/nThreads)+1;
	int tagID = currentWorker%nThreads;

	//MPI_Send(graphStr, cnt, MPI_CHAR, destination, 0, MPI_COMM_WORLD);
	MPI_Send(cstr, cnt, MPI_CHAR, rank, tagID, MPI_COMM_WORLD);
	if(Settings::debugMSG)
		cout<<"Master sends: "<<cstr<<" to worker: "<<currentWorker<<" ["<<rank<<","<<tagID<<"]"<<endl;
}

void Master::start(string fileName, bool ignoreFirstLine, int numMachine, int nThreads, bool approximate, int numPartitions, int nCorePoints, int minClusterSize) {

	this->nThreads = nThreads;
	int numWorkers = (numMachine-1)*nThreads;

	//load file from disk
	std::vector<Object*> allObjects = CVReader::read(fileName, ignoreFirstLine);

	int partitioningType = 0;//0- random partitions
							//1- smarter partitions to allow close objects partitioned together
	if(approximate)
		partitioningType = 1;

	//send its name to workers to load it there
	for(int i=1;i<numMachine;i++)
	{
		//Data Loading Task
		char fileName_[500];
		sprintf(fileName_, "f%s,%d,%d,%d,%d,%d,%d\n", fileName.c_str(), ignoreFirstLine, partitioningType,nCorePoints,Settings::debugMSG, numPartitions, nThreads);
		int cnt=strlen(fileName_)+1;
		MPI_Send(fileName_, cnt, MPI_CHAR, i, 0, MPI_COMM_WORLD);

		if(Settings::debugMSG)
			cout<<"Filename sent to worker: "<<i<<endl;
	}

	//waiting for all workers to load the data
	int tempNumWorkers = numMachine-1;
	while(tempNumWorkers>0)
	{
		MPI_Status status;
		char str_message[2000];
		MPI_Recv(str_message, 2000, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int srcThread = status.MPI_SOURCE*nThreads+status.MPI_TAG;

		if(Settings::debugMSG)
		{
			cout<<"Master received masg: "<<str_message<<endl;
		}

		if(str_message[0]=='l')
		{
			if(Settings::debugMSG)
				printf("Master received ack from %d(%d:%d)\n", status.MPI_SOURCE*nThreads+status.MPI_TAG, status.MPI_SOURCE, status.MPI_TAG);
		}

		tempNumWorkers--;
		if(Settings::debugMSG)
		{
			cout<<"Remaining workers to complete loading: "<<tempNumWorkers<<endl;
		}
	}

	//do the partitioning
	if(partitioningType==0) {

		if(Settings::debugMSG)
			cout<<"Partitioning data ... "<<flush;

		//generate random blocks
		rbs.createBlocks(numPartitions);
		//iterate over elements, assign random (equal size) elements to blocks
		rbs.assignElements(allObjects);
		rbs.prepareBlocks();
		rbs.prepareAllElements();

		if(Settings::debugMSG)
			cout<<"DONE."<<endl<<flush;
	}
	if(partitioningType==1) {

	}

	//send tasks to workers
	vector<std::pair<int, int>* > tasks;
	for(vector<RandomBlock*>::iterator iter1 = rbs.getBegin();iter1!=rbs.getEnd();iter1++) {

		//send single block
		tasks.push_back(new std::pair<int, int>((*iter1)->getID(), -1));

		//send pairs of blocks
		vector<RandomBlock*>::iterator iter2 = iter1;
		iter2++;
		for(;iter2!=rbs.getEnd();iter2++) {
			tasks.push_back(new std::pair<int, int>((*iter1)->getID(), (*iter2)->getID()));
		}
	}

	/*//number of pairs per worker
	int numTasksPerWorker = ceil(((double)tasks.size())/numWorkers);

	int currentWorker = nThreads;
	vector<std::pair<int, int>*> toSend;
	int numSentTasks = 0;
	for(vector<std::pair<int, int>* >::iterator iter = tasks.begin();iter!=tasks.end();iter++) {
		toSend.push_back(*iter);

		if(toSend.size()==numTasksPerWorker) {
			//send collected tasks to a worker
			std::stringstream msgStr;

			msgStr<<"t";
			for(int i=0;i<toSend.size();i++)
				msgStr<<","<<toSend[i]->first<<"|"<<toSend[i]->second;
			msgStr<<"\n";

			const std::string& tmp = msgStr.str();
			const char* cstr = tmp.c_str();

			int cnt=strlen(cstr);//+1;
			MPI_Status status;
			int rank = currentWorker/nThreads;
			int tagID = currentWorker%nThreads;

			//MPI_Send(graphStr, cnt, MPI_CHAR, destination, 0, MPI_COMM_WORLD);
			MPI_Send(cstr, cnt, MPI_CHAR, rank, tagID, MPI_COMM_WORLD);
			if(Settings::debugMSG)
				cout<<"Master sends: "<<cstr<<" to worker: "<<tagID<<endl;

			numSentTasks+=toSend.size();
			numTasksPerWorker = ceil(((double)tasks.size()-numSentTasks)/(numWorkers-currentWorker));
			currentWorker++;
			toSend.clear();
		}
	}*/

	bool* availableWorker = new bool[numWorkers];
	for(int i=0;i<numWorkers;i++)
	{
		availableWorker[i] = false;
		sendTask(tasks.back(), i);
		tasks.pop_back();
	}

	//Receive further messages
	int maxMsgSize = 10000000;
	char* str_message = new char[maxMsgSize];

	priority_queue<orderedLists> allOrderedLists;
	vector<orderedLists*> tempList;//a list used to collect and delete the created lists
	vector<list<PossibleEdge*>*> tempList2;//a list used to collect and delete the created lists
	while(true)
	{
		//wait for results from workers, it seems that I will need to wait for all results before building the final MST
		MPI_Status status;
		MPI_Recv(str_message, maxMsgSize, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int srcThread = (status.MPI_SOURCE-1)*nThreads+status.MPI_TAG;
		availableWorker[srcThread] = true;
		if(Settings::debugMSG)
			cout<<"Received results from worker "<<srcThread<<" with length = "<<strlen(str_message)<<endl;

		char command = str_message[0];
		switch(command) {
		case 'r':
			//send a new task before parsing
			if(tasks.size()>0) {
				//pick a worker, if possible
				int pickedWorker = -1;
				for(int i=0;i<numWorkers;i++) {
					if(availableWorker[i])
					{
						pickedWorker = i;
						break;
					}
				}

				if(pickedWorker>-1) {
					availableWorker[pickedWorker] = false;
					sendTask(tasks.back(), pickedWorker);
					tasks.pop_back();
				}
			}

			//parse results
			list<PossibleEdge*>* rList = parseEdgesMSG(str_message);
			tempList2.push_back(rList);

			orderedLists* ol = new orderedLists(rList);
			if(ol->getSize()>0) {
				allOrderedLists.push(*ol);
				tempList.push_back(ol);
			}
			else
				delete ol;

			if(Settings::debugMSG)
				cout<<"Master populated a list of resulted local MST edges with #edges = "<<rList->size()<<endl;

			if(Settings::debugMSG)
				cout<<"#remaining tasks = "<<tasks.size()<<endl;
			break;
		}

		//if all workers are available, break
		bool breakFromTheloop = true;
		for(int i=0;i<numWorkers;i++)
		{
			if(!availableWorker[i])
			{
				breakFromTheloop = false;
				break;
			}
		}

		if(breakFromTheloop) {
			if(Settings::debugMSG) {
				cout<<"All tasks are done, breaking to build the clusters."<<endl;
			}
			break;
		}
	}

	delete[] str_message;

	//builld the global MST and do clustering using this tree
	list<PossibleEdge*>* total = new list<PossibleEdge*>();

	//order all lists using an efficient technique
	while(allOrderedLists.size()>0) {
		orderedLists tempOL = allOrderedLists.top();
		allOrderedLists.pop();

		total->push_back(tempOL.getCurrentEdge());

		if(tempOL.advance())
			allOrderedLists.push(tempOL);
	}

	for(int i=0;i<tempList.size();i++) {
		delete tempList[i];
	}

	for(int i=0;i<tempList2.size();i++) {
		delete tempList2[i];
	}

	GraphX* mst = new GraphX();

	rbs.buildMSTFromListOfEdges(rbs.allElements, mst, total, false);

	Clusters* clusters = rbs.doClustering(mst, total, minClusterSize);
        std::ofstream outF (output.c_str(), std::ofstream::out);
        clusters->printElemValues(outF);

	delete mst;

	for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) {
		delete *iter;
	}
	delete total;

//	clusters->printElemValues(cout);
	cout<<"#Clusters = "<<clusters->getSize()<<endl;
	cout<<"#ExactClusters = "<<clusters->getSize()<<endl;
}
