/*
 * BlocksSystem.cpp
 *
 *  Created on: 9 Oct 2017
 *      Author: eabdelha
 */

#include<set>
#include<map>
#include<vector>
#include<list>
#include<queue>
#include <time.h>
#include <boost/pending/disjoint_sets.hpp>
#include "RandomBlocksSystem.h"
#include "GraphX.h"
#include "Hierarchy.h"
#include "CondenseHierarchy.h"
#include "PossibleEdge.h"
#include "TimeUtility.h"
#include "Utils.h"

//long RandomBlocksSystem::time1 = 0;

struct less_than_key
{
    inline bool operator() (PossibleEdge*& struct1, PossibleEdge*& struct2)
    {
        return (struct1->dist < struct2->dist);
    }
};

void getSortedFromSorted(list<PossibleEdge*>* l1, list<PossibleEdge*>* l2, list<PossibleEdge*>* r) {
	list<PossibleEdge*>::iterator iter1 = l1->begin();
	list<PossibleEdge*>::iterator iter2 = l2->begin();

	while(true) {

		if(iter1==l1->end() && iter2==l2->end())
			break;

		double dist1 = -1;
		double dist2 = -1;
		if(iter1!=l1->end())
			dist1 = (*iter1)->dist;
		if(iter2!=l2->end())
			dist2 = (*iter2)->dist;

		if(dist1!=-1) {
			if(dist2==-1 || dist1<=dist2) {
				r->push_back((*iter1));
				iter1++;
			}
		}
		if(dist2!=-1) {
			if(dist1==-1 || dist2<=dist1) {
				r->push_back((*iter2));
				iter2++;
			}
		}
	}
}

/**
 * add l1 elements to l2 sorted
 */
void getSortedFromSorted(list<PossibleEdge*>* l1, list<PossibleEdge*>* l2) {
	list<PossibleEdge*>::iterator iter1 = l1->begin();
	list<PossibleEdge*>::iterator iter2 = l2->begin();

	while(true) {

		if(iter1==l1->end() || (iter1==l1->end() && iter2==l2->end()))
			break;

		double dist1 = -1;
		double dist2 = -1;
		if(iter1!=l1->end())
			dist1 = (*iter1)->dist;
		if(iter2!=l2->end())
			dist2 = (*iter2)->dist;

		bool inserted = false;

		if(dist1!=-1) {
			if(dist2==-1 || dist1<=dist2) {
				l2->insert(iter2, (*iter1));
				iter1++;
				//iter2++;
				inserted = true;
			}
		}
		if(dist2!=-1 && !inserted) {
			if(dist1==-1 || dist2<=dist1) {
				iter2++;
			}
		}
	}
}

/**
 * add l1 elements to l2 sorted
 */
void getSortedFromSorted(vector<PossibleEdge*>* l1, list<PossibleEdge*>* l2) {
	vector<PossibleEdge*>::iterator iter1 = l1->begin();
	list<PossibleEdge*>::iterator iter2 = l2->begin();

	while(true) {

		if(iter1==l1->end() || (iter1==l1->end() && iter2==l2->end()))
			break;

		double dist1 = -1;
		double dist2 = -1;
		if(iter1!=l1->end())
			dist1 = (*iter1)->dist;
		if(iter2!=l2->end())
			dist2 = (*iter2)->dist;

		bool inserted = false;

		if(dist1!=-1) {
			if(dist2==-1 || dist1<=dist2) {
				l2->insert(iter2, (*iter1));
				iter1++;
				//iter2++;
				inserted = true;
			}
		}
		if(dist2!=-1 && !inserted) {
			if(dist1==-1 || dist2<=dist1) {
				iter2++;
			}
		}
	}
}

void RandomBlocksSystem::buildMSTFromListOfEdges(vector<Object*> objectsList, GraphX* mst, list<PossibleEdge*>* total, bool needMap) {
	int size = objectsList.size();

	std::vector<int>  rank (size+1);
	std::vector<int>  parent (size+1);
	boost::disjoint_sets<int*,int*> ds(&rank[0], &parent[0]);

	map<int, int> originalToSequenceMap;
	int count = 0;
	for(vector<Object*>::iterator iter = objectsList.begin();iter!=objectsList.end();iter++) {
		ds.make_set(count);
		if(needMap)
			originalToSequenceMap.insert(std::pair<int, int>((*iter)->getID(), originalToSequenceMap.size()));
		mst->AddNode((*iter)->getID(), 1, (*iter));
		count++;
	}

	for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();) {
		//pick the edge with smallest distance
		PossibleEdge* pe = *iter;

		//try to add this edge if possible
		int srcOrigID = pe->src->getID();
		int destOrigID = pe->dest->getID();
		if(needMap) {
			srcOrigID = originalToSequenceMap.find(pe->src->getID())->second;
			destOrigID = originalToSequenceMap.find(pe->dest->getID())->second;
		}
		if(ds.find_set(srcOrigID)!=ds.find_set(destOrigID)) {
			mst->addEdge(pe->src->getID(), pe->dest->getID(), pe->dist);
			ds.union_set(srcOrigID, destOrigID);
			iter++;
		} else {
			delete *iter;
			iter = total->erase(iter);
		}

		//break if all nodes are included in mst
		if(mst->getNumOfEdges()>=(mst->getNumOfNodes()-1) && mst->isConnected()) {
			while(iter!=total->end()) { PossibleEdge* tempPe = *iter; iter = total->erase(iter); delete(tempPe); }
			break;
		}
	}

//	cout<<"#Edges = "<<mst->getNumOfEdges()<<endl;
//	cout<<*mst<<endl;
}

RandomBlock* RandomBlocksSystem::getBlock(int index) {
	RandomBlock* rb = blocks[index];
	if(rb->getID()!=index)
		cout<<"Problemaaaa 11831991291032"<<endl;
	return rb;
}

Object* RandomBlocksSystem::getObject(int index) {
	Object* obj = allObjects[index];
	if(obj->getID()!=index)
		cout<<"Problemaaaa 39490327235. index = "<<index<<", obj->getID() = "<<obj->getID()<<endl;
	return obj;
}

void RandomBlocksSystem::createBlocks(int numPartitions) {
	for(int i=0;i<numPartitions;i++) {
		RandomBlock* rb = new RandomBlock();
		rb->setID(i);
		blocks.push_back(rb);
	}
}

void RandomBlocksSystem::assignElements(vector<Object*> allObjects) {
	//this->allObjects(allObjects);//.insert(allObjects.begin(), allObjects.end());
	this->allObjects.assign(allObjects.begin(), allObjects.end());

	int count = 0;
	for(std::vector<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {

		//assign elements
		blocks[count%blocks.size()]->addElement(*iter);
		count++;

		//add all elements to each block
		for(int i=0;i<blocks.size();i++) {
			blocks[i]->addObject(*iter);
		}
	}
}

void RandomBlocksSystem::prepareBlocks() {
	for(int i=0;i<blocks.size();i++) {
		blocks[i]->prepare();
	}
}

void RandomBlocksSystem::prepareSortedElements4Blocks() {
	for(int i=0;i<blocks.size();i++) {
		blocks[i]->prepareSortedElements();
	}
}

//void RandomBlocksSystem::addBlock(RandomBlock* block) {
//	this->blocks.push_back(block);
//}

void RandomBlocksSystem::computeCoreDistances(int nCorePoints) {
	for(vector<RandomBlock*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		(*iter)->computeCoreDistances(nCorePoints);
	}
}

void RandomBlocksSystem::print(ostream& out) {
	for(vector<RandomBlock*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		(*iter)->print(out);
	}
}

void RandomBlocksSystem::prepareAllElements() {
	for(vector<RandomBlock* >::iterator blocksIter = blocks.begin();blocksIter!=blocks.end();blocksIter++)
		for(vector<Object*>::iterator iter = (*blocksIter)->getElementsIter();iter!=(*blocksIter)->getElementsEndIter();iter++)
			allElements.push_back((*iter));
}

Clusters* RandomBlocksSystem::cluster(int minClusterSize) {
	//prepare for clustering
	GraphX* mst = new GraphX(1,0);
	list<PossibleEdge*>* total = new list<PossibleEdge*>();

	for(vector<RandomBlock* >::iterator blocksIter = blocks.begin();blocksIter!=blocks.end();blocksIter++)
		for(vector<Object*>::iterator iter = (*blocksIter)->getElementsIter();iter!=(*blocksIter)->getElementsEndIter();iter++)
			allElements.push_back((*iter));

	buildMST1(mst, total);
	Clusters* clusters = doClustering(mst, total, minClusterSize);

	delete mst;

	for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) {
		delete *iter;
	}
	delete total;

	return clusters;
}

Clusters* RandomBlocksSystem::clusterApprox1(int minClusterSize, double sampleRatio) {
	//prepare for clustering
	GraphX* mst = new GraphX(1,0);
	list<PossibleEdge*>* total = new list<PossibleEdge*>();

	for(vector<RandomBlock* >::iterator blocksIter = blocks.begin();blocksIter!=blocks.end();blocksIter++)
		for(vector<Object*>::iterator iter = (*blocksIter)->getElementsIter();iter!=(*blocksIter)->getElementsEndIter();iter++) {
			allElements.push_back((*iter));
		}

	buildMSTApprox1(mst, total, sampleRatio);
	Clusters* clusters = doClustering(mst, total, minClusterSize);

	for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) {
		delete *iter;
	}
	delete total;
	delete mst;

	return clusters;
}

void RandomBlocksSystem::buildMSTFromListofPairs(GraphX* mst, list<PossibleEdge*>* total, vector<std::pair<int, int>*> blocksPairs) {

	priority_queue<orderedLists> allOrderedLists;

	vector<orderedLists*> tempList;//a list used to collect and delete the created lists
	vector<list<PossibleEdge*>*> tempList2;//a list used to collect and delete the created lists
	for(int i=0;i<blocksPairs.size();i++) {
		std::pair<int, int>* currentTask = blocksPairs[i];
		list<PossibleEdge*>* l = buildLocalMSTofPairs(getBlock(currentTask->first), getBlock(currentTask->second));
		tempList2.push_back(l);

		orderedLists* ol = new orderedLists(l);
		if(ol->getSize()>0) {
			allOrderedLists.push(*ol);
			tempList.push_back(ol);
		}
		else
			delete ol;
	}
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
	//build the global MST
	buildMSTFromListOfEdges(allElements, mst, total, false);
}

void RandomBlocksSystem::buildMST1(GraphX* mst, list<PossibleEdge*>* total) {

	priority_queue<orderedLists> allOrderedLists;

	vector<orderedLists*> tempList;//a list used to collect and delete the created lists
	vector<list<PossibleEdge*>*> tempList2;//a list used to collect and delete the created lists

	int count = 0;
	for(vector<RandomBlock* >::iterator iter1 = blocks.begin();iter1!=blocks.end();iter1++) {
		//build local MST of pairs
		for(vector<RandomBlock* >::iterator iter2 = iter1+1;iter2!=blocks.end();iter2++) {

			list<PossibleEdge*>* l = buildLocalMSTofPairs((*iter1), (*iter2));
			tempList2.push_back(l);

			orderedLists* ol = new orderedLists(l);
			if(ol->getSize()>0) {
				allOrderedLists.push(*ol);
				tempList.push_back(ol);
			}
			else
				delete ol;
			//getSortedFromSorted(l, total);
			count++;
			cout<<count<<"/"<<(blocks.size()*(blocks.size()-1)/2)<<endl<<flush;
		}
	}

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

//	for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++)
//		cout<<(*iter)->dist<<endl;

	//build the global MST
	buildMSTFromListOfEdges(allElements, mst, total, false);
}

void RandomBlocksSystem::buildMSTApprox1(GraphX* mst, list<PossibleEdge*>* total, double sampleRatio) {

	srand (time(NULL));
	int count = 0;
	for(vector<RandomBlock* >::iterator iter1 = blocks.begin();iter1!=blocks.end();iter1++) {
		//build local MST of pairs
		for(vector<RandomBlock* >::iterator iter2 = iter1+1;iter2!=blocks.end();iter2++) {

			double r = ((double) rand() / (RAND_MAX));
			if(r<sampleRatio) {
				list<PossibleEdge*>* l = buildLocalMSTofPairs((*iter1), (*iter2));
				getSortedFromSorted(l, total);
				delete l;
				count++;
			}
		}
	}

	cout<<endl<<"The number of selected pairs = "<<count<<", total = "<<((blocks.size()*(blocks.size()-1))/2.0)<<endl;

	//build the global MST
	buildMSTFromListOfEdges(allElements, mst, total, false);
}

/**
 * complete the clustering given the MST graph
 */
Clusters* RandomBlocksSystem::doClustering(GraphX* mst, list<PossibleEdge*>* edgesSorted, int minClusterSize) {
//	cout<<"Clustering using the MST ... ";
	//create the hierarchy tree
	//root node of the hierarchy tree
	H_Node* root = NULL;
	double maxDistLink = 0;

	//map for the hierarchy tree
	map<int, H_Node*> hierarchy_map;//map to be used while creating
	std::vector<int>  rank1 (Object::getMaxID()+1);
	std::vector<int>  parent1 (Object::getMaxID()+1);
	boost::disjoint_sets<int*,int*> ds1(&rank1[0], &parent1[0]);
	for(vector<Object*>::iterator iter = allElements.begin();iter!=allElements.end();iter++) {
		Object* obj = (*iter);
		ds1.make_set(obj->getID());
	}

	int count = 0;
	for(list<PossibleEdge*>::iterator iter1 = edgesSorted->begin();iter1!=edgesSorted->end();iter1++)
	{
		count++;//cout<<count<<endl<<flush;
		Object* src = (*iter1)->src;
		NodeX* srcNode = mst->getNodeWithID(src->getID());
		Object* dest = (*iter1)->dest;
		NodeX* destNode = mst->getNodeWithID(dest->getID());
		double dist = (*iter1)->dist;

		H_Node* hnSrc;
		map<int, H_Node*>::iterator iterHM = hierarchy_map.find(ds1.find_set(src->getID()));
		if(iterHM==hierarchy_map.end())
		{
			hnSrc = new H_Node();
			hnSrc->setNode(srcNode);
		} else {
			hnSrc = iterHM->second;
		}

		H_Node* hnDest;
		iterHM = hierarchy_map.find(ds1.find_set(dest->getID()));
		if(iterHM==hierarchy_map.end())
		{
			hnDest = new H_Node();
			hnDest->setNode(destNode);
		} else {
			hnDest = iterHM->second;
		}

		H_Node* parent = new H_Node();
		parent->addNode(hnSrc, dist);hnSrc->addNode(parent, dist);
		parent->addNode(hnDest, dist);hnDest->addNode(parent, dist);

		if(dist>maxDistLink) {
			root = parent;
			maxDistLink = dist;
		}

		ds1.union_set(src->getID(), dest->getID());

		hierarchy_map.erase(ds1.find_set(src->getID()));
		hierarchy_map.erase(ds1.find_set(dest->getID()));
		hierarchy_map.insert(std::pair<int,H_Node*>(ds1.find_set(src->getID()), parent));
	}
	//print the hierarchy tree
	/*root->print();
	cout<<"Condensed tree:"<<endl;
	root->condense(2);
	root->print();*/

	//root->print();
	root->setNumItems();
//	root->print();
//	cout<<"Condensed advanced tree ..."<<endl;//<--- you have somthing wrong in your logic here!
	CH_Node* chn = new CH_Node();
	chn->setStart(0);
	chn->setWidth(root->getNumItems());
	root->condense(minClusterSize, chn);
	chn->traverseAreas();
//	chn->print(0);

//	cout<<"Clusters ...."<<endl;
	Clusters* clusters = new Clusters();
	vector<CH_Node*> clustersTemp;
	chn->getClusters(&clustersTemp);
	for(vector<CH_Node*>::iterator iter = clustersTemp.begin();iter!=clustersTemp.end();iter++) {
//		(*iter)->print(0, true);
//		(*iter)->printAsList();
		Cluster* cluster = new Cluster();
		cluster->populateObjects((*iter)->getDataNode());
		clusters->addClustrer(cluster);
	}
	CH_Node::deleteCh_Node(chn);
	root->deleteDescendents();
	delete root;

//	cout<<"DONE."<<endl;
	return clusters;
}

list<PossibleEdge*>* RandomBlocksSystem::buildLocalMSTofPairs(RandomBlock* b1, RandomBlock* b2) {
//	cout<<"Local MST for: "<<b1->id<<", "<<b2->id<<endl;

	list<PossibleEdge*>* total = new list<PossibleEdge*>();
	vector<PossibleEdge*>* l1 = b1->buildLocalMST_Prims();
	getSortedFromSorted(l1, total);
	delete l1;
//	cout<<"l1"<<endl;for(list<PossibleEdge*>::iterator iter = l1->begin();iter!=l1->end();iter++) (*iter)->print(cout);
//	cout<<"total"<<endl;for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) (*iter)->print(cout);
	vector<PossibleEdge*>* l2 = b2->buildLocalMST_Prims();
	getSortedFromSorted(l2, total);
	delete l2;
//	cout<<"l2"<<endl;for(list<PossibleEdge*>::iterator iter = l2->begin();iter!=l2->end();iter++) (*iter)->print(cout);
//	cout<<"total"<<endl;for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) (*iter)->print(cout);
	vector<PossibleEdge*>* l3 = buildBipartiteMSTofPairs_Fast(b1, b2);
	//list<PossibleEdge*>* l3 = buildBipartiteMSTofPairs(b1, b2);
	getSortedFromSorted(l3, total);
	delete l3;
//	cout<<"l3"<<endl;for(list<PossibleEdge*>::iterator iter = l3->begin();iter!=l3->end();iter++) (*iter)->print(cout);
//	cout<<"total"<<endl;for(list<PossibleEdge*>::iterator iter = total->begin();iter!=total->end();iter++) (*iter)->print(cout);
	vector<Object*> objectsList;
	for(vector<Object*>::iterator iter = b1->getElementsIter();iter!=b1->getElementsEndIter();iter++)	objectsList.push_back((*iter));
	for(vector<Object*>::iterator iter = b2->getElementsIter();iter!=b2->getElementsEndIter();iter++)	objectsList.push_back((*iter));

	GraphX* mst = new GraphX(1,0);
	buildMSTFromListOfEdges(objectsList, mst, total, true);
	delete mst;

	return total;
}

/**
 * build bipartite MST using the two partitions
 */
list<PossibleEdge*>* RandomBlocksSystem::buildBipartiteMSTofPairs(RandomBlock* b1, RandomBlock* b2) {
	//generate the list of possible edges
	map<int, PossibleEdge*> edges;
	for(vector<Object*>::iterator iter1 = b1->getElementsIter();iter1!=b1->getElementsEndIter();iter1++) {
		edges.insert(std::pair<int, PossibleEdge*>((*iter1)->getID(), new PossibleEdge((*iter1), NULL, -1)));
	}
	for(vector<Object*>::iterator iter2 = b2->getElementsIter();iter2!=b2->getElementsEndIter();iter2++) {
		edges.insert(std::pair<int, PossibleEdge*>((*iter2)->getID(), new PossibleEdge((*iter2), NULL, -1)));
	}

	//iterate over the first list, assign the minimum distances to both lists
	for(vector<Object*>::iterator iter1 = b1->getElementsIter();iter1!=b1->getElementsEndIter();iter1++) {
		Object* candidate1 = (*iter1);
		for(vector<Object*>::iterator iter2 = b2->getElementsIter();iter2!=b2->getElementsEndIter();iter2++) {
			Object* candidate2 = (*iter2);
			double dist = candidate1->coreDistanceTo(candidate2);

			PossibleEdge* peTemp;
			peTemp = edges.find(candidate1->getID())->second;
			if(peTemp->dist==-1 || dist<peTemp->dist) {
				peTemp->dist = dist;
				peTemp->dest = candidate2;
			}

			peTemp = edges.find(candidate2->getID())->second;
			if(peTemp->dist==-1 || dist<peTemp->dist) {
				peTemp->dist = dist;
				peTemp->dest = candidate1;
			}
		}
	}

	list<PossibleEdge*>* orderedEdges = new list<PossibleEdge*>();
	map<int, set<int>*> mapset;
//	GraphX mst = new GraphX();
	while(edges.size()>0) {
		double minDistance = edges.begin()->second->dist;
		PossibleEdge* minPE = edges.begin()->second;

		for(map<int, PossibleEdge*>::iterator iter = edges.begin();iter!=edges.end();iter++) {
			PossibleEdge* tempPE = iter->second;
			if(tempPE->dist<minDistance) {
				minDistance = tempPE->dist;
				minPE = tempPE;
			}
		}

//		mst.AddNode(minPE->src->getID(), 1);
//		mst.AddNode(minPE->dest->getID(), 1);
//		mst.addEdge(minPE->src->getID(), minPE->dest->getID(), minDistance);

		edges.erase(minPE->src->getID());

		int key1 = minPE->src->getID();
		int key2 = minPE->dest->getID();
		if(key1>key2) {
			int temp = key1;
			key1 = key2;
			key2 = temp;
		}

		map<int, set<int>*>::iterator tempIter = mapset.find(key1);
		if(tempIter==mapset.end()) {
			set<int>* temp = new set<int>();
			temp->insert(key2);
			mapset.insert(std::pair<int, set<int>*>(key1, temp));
			orderedEdges->push_back(minPE);
		} else {
			set<int>* tempSet = tempIter->second;
			if(tempSet->find(key2)==tempSet->end()) {
				tempSet->insert(key2);
				orderedEdges->push_back(minPE);
			}
		}
	}

//	cout<<mst<<endl;
	return orderedEdges;
}

/**
 * build bipartite MST using the two partitions
 */
vector<PossibleEdge*>* RandomBlocksSystem::buildBipartiteMSTofPairs_Fast(RandomBlock* b1, RandomBlock* b2) {

	vector<PossibleEdge*>* orderedEdges = new vector<PossibleEdge*>();
	//orderedEdges->reserve(b1->getNumberofElements()+b2->getNumberofElements());

	//iterate over the first list, get the minimum distances to the other list
	for(vector<Object*>::iterator iter1 = b1->getElementsIter();iter1!=b1->getElementsEndIter();iter1++) {
		Object* candidate1 = (*iter1);

		std::pair<Object* , double>* closest = b2->getClosestElement_CoreDistance(candidate1);

		PossibleEdge* pe = new PossibleEdge(candidate1, closest->first, closest->second);
		delete closest;

		orderedEdges->push_back(pe);
	}

	//iterate over the second list, get the minimum distances to the other list
	for(vector<Object*>::iterator iter2 = b2->getElementsIter();iter2!=b2->getElementsEndIter();iter2++) {
		Object* candidate2 = (*iter2);

		std::pair<Object* , double>* closest = b1->getClosestElement_CoreDistance(candidate2);

		PossibleEdge* pe = new PossibleEdge(candidate2, closest->first, closest->second);
		delete closest;

		orderedEdges->push_back(pe);
	}

	std::sort(orderedEdges->begin(), orderedEdges->end(), less_than_key());

	return orderedEdges;
}

/*void RandomBlocksSystem::buildMST() {

	vector<Object*> allObjects;

	//collect the unique objects
	for(vector<RandomBlock*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		RandomBlock* temp = (*iter);
		vector<Object*> tempList = temp->getElements();
//		allObjects.insert(tempList->begin(), tempList->end());
		std::copy(tempList.begin(), tempList.end(), std::back_insert_iterator<std::vector<Object*> >(allObjects));
	}

	//create the disjoint set thing for the Kruskal algorithm
	//////////////////

	std::vector<int>  rank (Object::getMaxID()+1);
	std::vector<int>  parent (Object::getMaxID()+1);
	boost::disjoint_sets<int*,int*> ds(&rank[0], &parent[0]);

	////////////////////////////////////

	GraphX* mst = new GraphX(1,0);
	for(vector<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {
		Object* obj = (*iter);
		mst->AddNode(obj->getID(), 1, obj);

		// for dsets
		ds.make_set(obj->getID());
	}

	vector<PossibleEdge*> finalAddedEdgesSorted;

	//create the minimum spanning tree
	double min = 0;
	double range = 10000;//10000;
	while(true) {

		cout<<"min = "<<min<<", max = "<<(min+range)<<endl;

		//sort objects based on the distance ascendingly
		list<PossibleEdge*> allObjectsSorted;

		cout<<"Get distances ..."<<flush;

		int count = 0;
		for(vector<RandomBlock* >::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
			RandomBlock* rb = (*iter);
			cout<<"finding distances ..."<<flush;
			vector<PossibleEdge*> tempDistances = rb->getPairsWithDistBound(min, min+range);
			cout<<"DONE."<<endl<<flush;

			cout<<"Sorting objects ... all.size = "<<allObjectsSorted.size()<<", partial.size = "<<tempDistances.size()<<flush;

			//apply merge sort
			list<PossibleEdge*>::iterator allElementsIter = allObjectsSorted.begin();
			vector<PossibleEdge*>::iterator partialListIter = tempDistances.begin();
			while(true) {
				if(allElementsIter==allObjectsSorted.end()) {
					allObjectsSorted.insert(allElementsIter, partialListIter, tempDistances.end());
					break;
				}
				if(partialListIter==tempDistances.end()) {
					break;
				}

				if((*allElementsIter)->dist<(*partialListIter)->dist) {
					allElementsIter++;
				} else {
					allElementsIter = allObjectsSorted.insert(allElementsIter, (*partialListIter));
					allElementsIter++;
					partialListIter++;
				}
			}

			//put in place
//			for(vector<PossibleEdge*>::iterator iter1 = tempDistances.begin();iter1!=tempDistances.end();iter1++) {
//
//				PossibleEdge* pe = (*iter1);
//
//				vector<PossibleEdge*>::iterator innerIter = allObjectsSorted.begin();
//				for(;innerIter!=allObjectsSorted.end();innerIter++) {
//					if(pe->dist<(*innerIter)->dist)
//						break;
//				}
//
//				allObjectsSorted.insert(innerIter, pe);
//			}
			cout<<"DONE."<<endl<<flush;
			count++;
			cout<<count<<"/"<<blocks.size()<<endl;
		}

		cout<<"Number of new recieved edges = "<<": "<<allObjectsSorted.size()<<endl;
		cout<<"DONE."<<flush;

		int numPrevEdges = mst->getNumOfEdges();
		cout<<"Generate MST ..."<<flush;
		//generate the mst tree
		list<PossibleEdge*>::iterator innerIter = allObjectsSorted.begin();
		for(;innerIter!=allObjectsSorted.end();innerIter++) {
//			(*innerIter)->print(cout);

			Object* src = (*innerIter)->src;
			NodeX* srcNode = mst->getNodeWithID(src->getID());
			Object* dest = (*innerIter)->dest;
			NodeX* destNode = mst->getNodeWithID(dest->getID());
			double dist = (*innerIter)->dist;

			if(ds.find_set(src->getID())==ds.find_set(dest->getID()))
			{
				continue;
			}

			//add an edge to the minimum spanning tree
			mst->addEdge(srcNode->getID(), destNode->getID(), dist);
			ds.union_set(src->getID(), dest->getID());

			vector<PossibleEdge*>::iterator iter1 = finalAddedEdgesSorted.begin();
			for(;iter1!=finalAddedEdgesSorted.end();iter1++) {
				if((*iter1)->dist>(*innerIter)->dist)
					break;
			}
			finalAddedEdgesSorted.insert(iter1, (*innerIter));
		}
		int numNewEdges = mst->getNumOfEdges() - numPrevEdges;
		cout<<"num added new edges = "<<numNewEdges<<". ";
		cout<<"DONE."<<flush;

		//check connected components
		map<int, int> componentIDToCountMap;
		for(map<int,NodeX*>::const_iterator nodesIter = mst->getNodesIterator();nodesIter!=mst->getNodesEndIterator();nodesIter++) {
			int nodeID = nodesIter->first;

			int componentID = ds.find_set(nodeID);

			if(componentIDToCountMap.find(componentID)==componentIDToCountMap.end()) {
				componentIDToCountMap.insert(std::pair<int, int>(componentID, 1));
			} else {
				int currentCount = componentIDToCountMap.find(componentID)->second;
				componentIDToCountMap.erase(componentID);
				componentIDToCountMap.insert(std::pair<int, int>(componentID, currentCount+1));
			}
		}

		cout<<"Number of connected components: "<<componentIDToCountMap.size()<<endl;
		for(map<int, int>::iterator tempIter = componentIDToCountMap.begin();tempIter!=componentIDToCountMap.end();tempIter++) {
			cout<<"\t"<<tempIter->first<<", "<<tempIter->second<<endl;
		}

		//check if it is connected
		if(mst->isConnected()) {
//			cout<<*mst;
			break;
		}

//		cout<<"SHADOW SMT:"<<*mst_shadow<<endl;
//		cout<<"SMT:"<<*mst<<endl;

		min = min + range;
	}

	//create the hierarchy tree
	//root node of the hirarchy tree
	H_Node* root = NULL;
	double maxDistLink = 0;

	//map for the hierarchy tree
	map<int, H_Node*> hierarchy_map;//map to be used while creating
	std::vector<int>  rank1 (Object::getMaxID()+1);
	std::vector<int>  parent1 (Object::getMaxID()+1);
	boost::disjoint_sets<int*,int*> ds1(&rank1[0], &parent1[0]);

	for(vector<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {
		Object* obj = (*iter);
		ds1.make_set(obj->getID());
	}

	for(vector<PossibleEdge*>::iterator iter1 = finalAddedEdgesSorted.begin();iter1!=finalAddedEdgesSorted.end();iter1++)
	{
		Object* src = (*iter1)->src;
		NodeX* srcNode = mst->getNodeWithID(src->getID());
		Object* dest = (*iter1)->dest;
		NodeX* destNode = mst->getNodeWithID(dest->getID());
		double dist = (*iter1)->dist;

		H_Node* hnSrc;
		map<int, H_Node*>::iterator iterHM = hierarchy_map.find(ds1.find_set(src->getID()));
		if(iterHM==hierarchy_map.end())
		{
			hnSrc = new H_Node();
			hnSrc->setNode(srcNode);
		} else {
			hnSrc = iterHM->second;
		}

		H_Node* hnDest;
		iterHM = hierarchy_map.find(ds1.find_set(dest->getID()));
		if(iterHM==hierarchy_map.end())
		{
			hnDest = new H_Node();
			hnDest->setNode(destNode);
		} else {
			hnDest = iterHM->second;
		}

		H_Node* parent = new H_Node();
		parent->addNode(hnSrc, dist);hnSrc->addNode(parent, dist);
		parent->addNode(hnDest, dist);hnDest->addNode(parent, dist);

//		cout<<"*** "<<parent<<","<<hnSrc<<","<<hnDest<<endl;

		if(dist>maxDistLink) {
			root = parent;
			maxDistLink = dist;
		}

		ds1.union_set(src->getID(), dest->getID());

		hierarchy_map.erase(ds1.find_set(src->getID()));
		hierarchy_map.erase(ds1.find_set(dest->getID()));
		hierarchy_map.insert(std::pair<int,H_Node*>(ds1.find_set(src->getID()), parent));
	}

	//print the hierarchy tree
//	root->print();
//	cout<<"Condensed tree:"<<endl;
//	root->condense(2);
//	root->print();

	//root->print();

	root->setNumItems();
//	root->print();
}*/

RandomBlocksSystem::~RandomBlocksSystem() {
	cout<<"calling RBS destructor ..."<<endl;

	for(int i=0;i<blocks.size();i++) {
		delete blocks[i];
	}

	cout<<"desctructor done!"<<endl;
}
