/*
 * Worker.cpp
 *
 *  Created on: 7 Feb 2018
 *      Author: eabdelha
 */

#include <mpi.h>
#include <boost/thread.hpp>
#include "CSVReader.h"
#include "Worker.h"
#include "Settings.h"
#include "Utils.h"
#include "PossibleEdge.h"

using namespace std;

std::pair<int, int>* Worker::parseTaskMSG(char* str_message) {
	pair<int, int>* task;
	int pos = 1;
	//get the first number of th epair
	int pos1 = getPos(str_message, '|', pos+1);

	char block1Str[10];
	strncpy (block1Str, str_message+pos+1, pos1-pos);
	int block1 = atoi(block1Str);

	pos = getPos(str_message, ',', pos1+1);

	char block2Str[10];
	if(pos==-1)
		strncpy (block2Str, str_message+pos1+1, strlen(str_message)-pos1);
	else
		strncpy (block2Str, str_message+pos1+1, pos-pos1);
	int block2 = atoi(block2Str);

	task = new pair<int, int>(block1, block2);

	return task;
}

void Worker::start(int rank) {
	this->machineID = rank;
	//this->nThreads = nThreads;
	char str_message[BUFFERSIZE];
	boost::thread** threads;

	//initialize
	//get filename from master
	if(Settings::debugMSG)
		cout<<"Worker is waiting for messages ..."<<endl;

	MPI_Status status;
	MPI_Recv(str_message, BUFFERSIZE, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	if(Settings::debugMSG)
		cout<<"Worker received load msg:"<<str_message<<endl<<flush;

	if(str_message[0]=='f')
	{
		//parse file name
		int pos = getPos(str_message, ',');
		char fileName[500];
		strncpy (fileName, str_message+1, pos-1);
		fileName[pos-1] = '\0';

		//parse the ignoreFirstLine bool
		int pos1 = getPos(str_message, ',', pos+1);
		char ignoreFLStr[10];
		strncpy (ignoreFLStr, str_message+pos+1, pos1-pos);
		int ignoreFirstLine = atoi(ignoreFLStr);

		//parse the partition type
		pos = getPos(str_message, ',', pos1+1);
		char partTypeStr[10];
		strncpy (partTypeStr, str_message+pos1+1, pos-pos1);
		int partitionType = atoi(partTypeStr);

		//parse number of core points
		pos1 = getPos(str_message, ',', pos+1);
		char nCoreStr[10];
		strncpy (nCoreStr, str_message+pos+1, pos1-pos);
		nCorePoints = atoi(nCoreStr);

		//parse the show debug messages flag
		pos = getPos(str_message, ',', pos1+1);
		char showDMStr[10];
		strncpy (showDMStr, str_message+pos1+1, pos-pos1);
		int showDM = atoi(showDMStr);
		Settings::debugMSG = false;
		if(showDM==1)
			Settings::debugMSG = true;

		//parse the number of partitions
		pos1 = getPos(str_message, ',', pos+1);
		char numPartsStr[10];
		strncpy (numPartsStr, str_message+pos+1, pos1-pos);
		int numPartitions = atoi(numPartsStr);

		//parse the number of threads
		char numThreadsStr[10];
		strncpy (numThreadsStr, str_message+pos1+1, strlen(str_message)-pos1);
		this->nThreads = atoi(numThreadsStr);

		if(Settings::debugMSG)
			printf("Worker %d: Received File Name: %s, ignore first line = %d, nCorePoints = %d, partition Type = %d, number of partitions = %d, number of threads = %d \n", machineID, fileName, ignoreFirstLine, nCorePoints, partitionType, numPartitions, nThreads);

		std::vector<Object*> allObjects = CVReader::read(fileName, ignoreFirstLine==1);

		if(Settings::debugMSG)
			printf("Data loaded at worker: %d\n",machineID);

		//do the partitioning
		if(partitionType==0) {

			if(Settings::debugMSG)
				cout<<"Partitioning data [Worker] ... "<<flush;

			//generate random blocks
			rbs.createBlocks(numPartitions);
			//iterate over elements, assign random (equal size) elements to blocks
			rbs.assignElements(allObjects);
			rbs.prepareBlocks();

			rbs.computeCoreDistances(nCorePoints);

			if(Settings::debugMSG)
				cout<<"DONE [Worker]."<<endl<<flush;
		}
		if(partitionType==1) {

		}

		rbs.prepareAllElements();

		//send notification to master
		char str_message_1[BUFFERSIZE];

		sprintf(str_message_1, "l");
		int cnt=strlen(str_message_1)+1;

		//initialize threads
		threads = new boost::thread*[nThreads];
		for(int i=0;i<nThreads;i++)
		{
			try
			{
				threads[i] = new boost::thread(&Worker::processing, this, i, &rbs);
			}
			catch(boost::thread_exception ex)
			{
				cout<<ex.what()<<endl;
				cout<<"Thread Exception"<<endl<<flush;
			}
		}

		MPI_Send(str_message_1, cnt, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	else
	{
		cout<<"uqskgkjsdfk"<<endl;
		exit(0);
	}

	for(int i=0;i<nThreads;i++)
	{
		threads[i]->join();
	}
}

void *Worker::processing(Worker* worker, int id, RandomBlocksSystem* rbs)
{
	usleep(1000000);
	int threadID = (worker->machineID-1)*worker->nThreads+id;

	char* str_message = new char[BUFFERSIZE+1];

	while(1)
	{
		//receive task
		if(Settings::debugMSG)
			cout<<"Worker "<<threadID<<" is waiting for messages ..."<<endl;

		MPI_Status status;
		MPI_Recv(str_message, BUFFERSIZE, MPI_CHAR, 0, id, MPI_COMM_WORLD, &status);

		if(Settings::debugMSG)
			cout<<"Worker received msg:"<<str_message<<endl<<flush;

		char command = str_message[0];
		switch(command) {
		case 't':
			std::pair<int, int>* task = parseTaskMSG(str_message);

			vector<PossibleEdge*>* l3;
			if(task->second!=-1) {
				RandomBlock* b1 = rbs->getBlock(task->first);
				RandomBlock* b2 = rbs->getBlock(task->second);
				l3 = rbs->buildBipartiteMSTofPairs_Fast(b1, b2);
			} else {
				RandomBlock* b1 = rbs->getBlock(task->first);
				l3 = b1->buildLocalMST_Prims();
			}

			delete task;

			//send results back to the master
			std::stringstream msgStr;
			msgStr<<"r";
			for(vector<PossibleEdge*>::iterator iter = l3->begin();iter!=l3->end();iter++) {
				//send collected tasks to the master
				msgStr<<","<<(*iter)->src->getID()<<"|"<<(*iter)->dest->getID()<<"|"<<(*iter)->dist;
			}
			msgStr<<'\0';

			const std::string& tmp = msgStr.str();
			const char* cstr = tmp.c_str();

			int cnt=strlen(cstr)+1;

			MPI_Send(cstr, cnt, MPI_CHAR, 0, id, MPI_COMM_WORLD);
			if(Settings::debugMSG)
				cout<<"Worker "<<threadID<<" sends results back to master. size = "<<cnt<<endl;

			for(int i=0;i<l3->size();i++)
				delete (*l3)[i];
			delete l3;

			break;
		}
	}

	delete[] str_message;

	while(1) {

	}
}
