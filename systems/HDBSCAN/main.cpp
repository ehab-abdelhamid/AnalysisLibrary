/*
 * main.cpp
 *
 *  Created on: 27 Sep 2017
 *      Author: eabdelha
 */

#include <fstream>
#include<iostream>
#include<cstdlib>
#include<mpi.h>
#include "CSVReader.h"
#include "Object.h"
#include "Utils.h"
#include "Block.h"
#include "BlocksSystem.h"
//#include "ClosePartition.h"
#include "RandomBlocksSystem.h"
#include "RandomBlock.h"
#include "TimeUtility.h"
#include "Master.h"
#include "Worker.h"
#include "Settings.h"

using namespace std;

long RandomBlock::time1 = 0;

bool Settings::debugMSG;

int Worker::nCorePoints = -1;

int main( int argc, char *argv[])
{
	//testing vector on line projection
	if(false) {
		Object* x = new Object();
		x->addValue(6); x->addValue(5);

		Object* l1 = new Object();
		l1->addValue(1); l1->addValue(1);

		Object* l2 = new Object();
		l2->addValue(5); l2->addValue(2);

		Block b;
		Object* r = b.getProjection(x, new pair<Object*, Object*>(l1, l2));
		r->print(cout);

		return 0;
	}

	if(false) {
		//1:(2,2), 2:(1.5,2), 3:(2,1.5), 4:(1.75,1.75), 5:(3,2), 6:(3,3), 7:(3,1), 8:(3,1.1)
		Object* o1 = new Object();o1->addValue(2);o1->addValue(2);
		Object* o2 = new Object();o2->addValue(1.9);o2->addValue(2);
		Object* o3 = new Object();o3->addValue(2);o3->addValue(1.9);
		Object* o4 = new Object();o4->addValue(1.5);o4->addValue(2);
		Object* o5 = new Object();o5->addValue(3);o5->addValue(2);
		Object* o6 = new Object();o6->addValue(3);o6->addValue(3);
		Object* o7 = new Object();o7->addValue(3);o7->addValue(1.2);
		Object* o8 = new Object();o8->addValue(3);o8->addValue(1.1);
		Block* b1 = new Block();
		b1->addElement(o1, true);b1->addElement(o2, true);b1->addElement(o3, true);b1->addElement(o4, true);
		b1->addElement(o5, true);b1->addElement(o6, true);b1->addElement(o7, true);b1->addElement(o8, true);

		BlocksSystem bs;
		bs.addBlock(b1);

		/*Object* o1 = new Object();o1->addValue(1);o1->addValue(1);
		Object* o2 = new Object();o2->addValue(1);o2->addValue(2);
		Object* o3 = new Object();o3->addValue(2);o3->addValue(3);
		Object* o4 = new Object();o4->addValue(4);o4->addValue(1);
		Object* o5 = new Object();o5->addValue(4);o5->addValue(3);
		Object* o6 = new Object();o6->addValue(6);o6->addValue(3);
		Object* o7 = new Object();o7->addValue(7);o7->addValue(1);
		Object* o8 = new Object();o8->addValue(7);o8->addValue(2);
		Object* o9 = new Object();o9->addValue(7);o9->addValue(4);
		Object* o10 = new Object();o10->addValue(8);o10->addValue(1);
		Object* o11 = new Object();o11->addValue(8);o11->addValue(2);
		Object* o12 = new Object();o12->addValue(9);o12->addValue(5);
		Object* o13 = new Object();o13->addValue(10);o13->addValue(2);

		Block* b1 = new Block();
		b1->addElement(o1);b1->addElement(o2);b1->addElement(o3);b1->addElement(o4);
		Block* b2 = new Block();b2->addElement(o3);b2->addElement(o4);b2->addElement(o5);b2->addElement(o6);b2->addElement(o7);b2->addElement(o8);b2->addElement(o9);
		Block* b3 = new Block();b3->addElement(o6);b3->addElement(o7);b3->addElement(o8);b3->addElement(o9);b3->addElement(o10);b3->addElement(o11);b3->addElement(o12);b3->addElement(o13);

		BlocksSystem bs;
		bs.addBlock(b1);bs.addBlock(b2);bs.addBlock(b3);*/

//		pair<Object*, double>* p1 = o1->getKShortestRechability(NULL, &bs, 3, 1);
//		cout<<p1->first->getID()<<", "<<p1->second<<endl;
//		pair<Object*, double>* p2 = o2->getKShortestRechability(NULL, &bs, 3, 1);
//		cout<<p2->first->getID()<<", "<<p2->second<<endl;

		bs.buildMST();

//		int coreK = 5;
//		int k = 12;
//		pair<Object*,double>* p = o5->getKShortestRechability(NULL, &bs, coreK, k);
//		cout<<"Here is the "<<k<<"th shortest reachability object: "<<endl;
//		p->first->print(cout);
//		cout<<"With diatnce: "<<p->second<<endl;

		cout<<"Done!"<<endl;

		return 0;
	}

	if(false) {
		//1:(2,2), 2:(1.5,2), 3:(2,1.5), 4:(1.75,1.75), 5:(3,2), 6:(3,3), 7:(3,1), 8:(3,1.1)
		Object* o1 = new Object();o1->addValue(2);o1->addValue(2);
		Object* o2 = new Object();o2->addValue(1.5);o2->addValue(2);
		Object* o3 = new Object();o3->addValue(2);o3->addValue(1.5);
		Object* o4 = new Object();o4->addValue(1.75);o4->addValue(1.75);
		Object* o5 = new Object();o5->addValue(3);o5->addValue(2);
		Object* o6 = new Object();o6->addValue(3);o6->addValue(3);
		Object* o7 = new Object();o7->addValue(3);o7->addValue(1);
		Object* o8 = new Object();o8->addValue(3);o8->addValue(1.1);
		RandomBlock* b1 = new RandomBlock();
		b1->addElement(o1);b1->addElement(o2);b1->addElement(o3);b1->addElement(o4);
		b1->addElement(o5);b1->addElement(o6);b1->addElement(o7);b1->addElement(o8);

		b1->buildLocalMST_Prims();
		cout<<"done!"<<endl;
		return 1;
	}

	//latest 4-1-2018
	if(false) {
		//O -- 1:(1,2), 2:(2,5), 3:(2.5,4), 4:(2,2), 5:(2,1)
		Object* o1 = new Object();o1->addValue(1);o1->addValue(2);
		Object* o2 = new Object();o2->addValue(2);o2->addValue(5);
		Object* o3 = new Object();o3->addValue(2.5);o3->addValue(4);
		Object* o4 = new Object();o4->addValue(2);o4->addValue(2);
		Object* o5 = new Object();o5->addValue(2);o5->addValue(1);
		RandomBlock* b1 = new RandomBlock();b1->setID(1);
		b1->addElement(o1);b1->addElement(o2);b1->addElement(o3);b1->addElement(o4);b1->addElement(o5);

		//X -- 1:(4,1), 2:(4,5), 3:(2,2.5), 4:(3,4), 5:(2,6), 6:(1.5,2)
		Object* x1 = new Object();x1->addValue(4);x1->addValue(1);
		Object* x2 = new Object();x2->addValue(4);x2->addValue(5);
		Object* x3 = new Object();x3->addValue(2);x3->addValue(2.5);
		Object* x4 = new Object();x4->addValue(3);x4->addValue(4);
		Object* x5 = new Object();x5->addValue(2);x5->addValue(6);
		Object* x6 = new Object();x6->addValue(1.5);x6->addValue(2);
		RandomBlock* b2 = new RandomBlock();b2->setID(2);
		b2->addElement(x1);b2->addElement(x2);b2->addElement(x3);b2->addElement(x4);b2->addElement(x5);b2->addElement(x6);

		RandomBlocksSystem* rbs = new RandomBlocksSystem();
		rbs->addBlock(b1);
		rbs->addBlock(b2);

		rbs->cluster(2);

		cout<<"done!"<<endl;
		return 1;
	}

	bool parallel = true;
	int numPartitions = 10;
//	string fileName = "../../Research/Code/Results/Dataset/Brain/ADNIMERGE_2.csv"; bool ignoreFirstLine = true;
	string fileName = "/data/Research/Code/Results/Dataset/Higgs/Higgs-training-nolabels.csv"; bool ignoreFirstLine = false;
//	string fileName = "/data/Research/Code/Results/Dataset/Higgs/Higgs-testing-nolabels.csv"; bool ignoreFirstLine = false;
//	string fileName = "/data/Research/Code/Results/Dataset/Credit/credit-card.csv"; bool ignoreFirstLine = true;
//	string fileName = "/data/Research/Code/out.csv"; bool ignoreFirstLine = false;
//	string fileName = "/data/Research/Code/out-1D.csv"; bool ignoreFirstLine = false;

	int nCorePoints = 5;
	int minClusterSize = 5;
	bool exact = true;//these are connected some how!
	bool approximate = false;

	//get the show debug messages flag
	char * argShowDebugMsg = getCmdOption(argv, argv + argc, "-showDebugMSG");
	if(argShowDebugMsg)
	{
		if(argShowDebugMsg[0]=='1')
			Settings::debugMSG = true;
		else
			Settings::debugMSG = false;
	}
	else
		Settings::debugMSG = false;
	//get the number of threads
	char * argPartitions = getCmdOption(argv, argv + argc, "-nPartitions");
	if(argPartitions)
	{
		numPartitions = atoi(argPartitions);
	}
	//get the filename
	char * argFilename = getCmdOption(argv, argv + argc, "-file");
	if(argFilename)
	{
		fileName = argFilename;
	}
        //get the output filename
        string outFileName;
        char * argOutFilename = getCmdOption(argv, argv + argc, "-output");
        if(argOutFilename)
        {
                outFileName = argOutFilename;
        }


	//2. check why your results against the Python system
	//3. check how to run your system on a distrubuted environment

	if(parallel) {
		//Initialise MPI
		int i;
		int rank;
		int numWorkers;
		int nThreads = 1;
		MPI_Status    status;
		srand (time(NULL));

		int provided = MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
		if(provided != MPI_THREAD_MULTIPLE){
			cout<<"Provided = "<<provided<<endl;
			cout<<"ERROR: Cannot initialise MPI with the desire level of support";
		}
		rank = MPI::COMM_WORLD.Get_rank();
		MPI_Comm_size(MPI_COMM_WORLD, &numWorkers);

		//load data from the file. Note that data is replicated.

		int cnt=0;
		if(rank==0)
		{
			//the master process
			cout<<"Master started ..."<<endl;

			//get the number of threads
			char * argThreads = getCmdOption(argv, argv + argc, "-threads");
			if(argThreads)
			{
				nThreads = atoi(argThreads);
			}

			//start the master node
			TimeUtility::StartCounterMicro();

			Master m;
                        m.output = outFileName;
			m.start(fileName, ignoreFirstLine, numWorkers, nThreads, approximate, numPartitions, nCorePoints, minClusterSize);

			double elapsed = TimeUtility::GetCounterMicro();
			cout<<"[EXACT] elapsed time in microseconds = "<<elapsed<<endl<<flush;
		}
		else
		{
			//start workers
			cout<<"Worker: "<<rank<<" ..."<<endl;

			Worker w;
			w.start(rank);
		}
	} else {
		std::vector<Object*> allObjects = CVReader::read(fileName, ignoreFirstLine);

		cout<<"Partitioning data ... "<<flush;

		RandomBlocksSystem rbs;

		//generate random blocks
		rbs.createBlocks(numPartitions);

		//iterate over elements, assign random (equal size) elements to blocks
		rbs.assignElements(allObjects);

		rbs.prepareBlocks();

		cout<<"DONE."<<endl<<flush;

		Clusters* clusters = NULL;
		Clusters* approxClusters = NULL;

		if(exact) {
			cout<<"Exact solution starts ..."<<endl<<flush;
			TimeUtility::StartCounterMicro();
			rbs.computeCoreDistances(nCorePoints);
			double elapsed = TimeUtility::GetCounterMicro();
			cout<<"[EXACT CORE DISTANCES] elapsed time in microseconds = "<<elapsed<<endl<<flush;

			rbs.prepareSortedElements4Blocks();

			TimeUtility::StartCounterMicro();
			clusters = rbs.cluster(minClusterSize);
			clusters->printElemValues(cout);
			cout<<"#Clusters = "<<clusters->getSize()<<endl;
			cout<<"#ExactClusters = "<<clusters->getSize()<<endl;
			elapsed = TimeUtility::GetCounterMicro();
			cout<<"[EXACT] elapsed time in microseconds = "<<elapsed<<endl<<flush;
		}

		if(approximate) {
			//the second approximate clustering
			cout<<"Approximate solution starts ..."<<endl<<flush;
			TimeUtility::StartCounterMicro();
			rbs.computeCoreDistances(nCorePoints);//this is repeated here for comparison purposes only
			double elapsed = TimeUtility::GetCounterMicro();
			cout<<"[APPROX CORE DISTANCES] elapsed time in microseconds = "<<elapsed<<endl;

			TimeUtility::StartCounterMicro();
			approxClusters = rbs.clusterApprox1(minClusterSize, 0.1);
		//		approxClusters->printElemValues(cout);
			cout<<"#ApproxClusters = "<<approxClusters->getSize()<<endl;
			elapsed = TimeUtility::GetCounterMicro();
			cout<<"[APPROX] elapsed time in microseconds = "<<elapsed<<endl;
			if(exact)
				cout<<"Purity = "<<approxClusters->computePurity(clusters)<<endl;
			delete approxClusters;
		}

		if(exact)
			delete clusters;

		cout<<"time1 = "<<RandomBlock::time1<<endl;

		//remove objects
		for(int i=0;i<allObjects.size();i++)
			delete allObjects[i];

		annClose();

		cout<<"All DONE."<<endl;

		return 1;
	}

	/*std::vector<Object*> sampledObjects;
	sampledObjects = allObjects;//sample(allObjects, (float)0.1);
	cout<<"#sampled objects = "<<sampledObjects.size()<<endl;

	vector<Block*>* smallestBlocks = new vector<Block*>();
	Block* block = new Block();
	block->addElemets(sampledObjects, true);
	block->partition(smallestBlocks);

	BlocksSystem bs;

	//add to a blocks System
	for(vector<Block*>::iterator iter = smallestBlocks->begin();iter!=smallestBlocks->end();iter++) {
		(*iter)->print(cout);
		bs.addBlock((*iter));
	}

	bs.prepareBlocks();
	bs.print(cout);

//	if(true) return 1;

	bs.buildMST();


//	std::ifstream file("../../Research/Code/Results/Dataset/Brain/ADNIMERGE_2.csv");
//
//	vector<Object*> allObjects;
//	CSVRow row;
//	while(file >> row)
//	{
//		Object* o = new Object();
//		for(unsigned int i=0;i<row.size();i++) {
//			o->addValue((float)row[i]);
//		}
//		o->print(cout);
//		allObjects.push_back(o);
//	}
//
//	cout<<"#objects = "<<allObjects.size()<<endl;

	cout<<"DONE!";*/
}
