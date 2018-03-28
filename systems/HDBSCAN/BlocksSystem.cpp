/*
 * BlocksSystem.cpp
 *
 *  Created on: 9 Oct 2017
 *      Author: eabdelha
 */

#include<set>
#include<map>
#include <boost/pending/disjoint_sets.hpp>
#include "BlocksSystem.h"
#include "GraphX.h"
#include "Hierarchy.h"
#include "CondenseHierarchy.h"
#include "PossibleEdge.h"

void BlocksSystem::addBlock(Block* block) {
	this->blocks.push_back(block);
}

void BlocksSystem::computeCoreDistances() {
	for(vector<Block*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		(*iter)->computeCoreDistances(this);
	}
}

void BlocksSystem::print(ostream& out) {
	for(vector<Block*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		(*iter)->print(out, false);
	}
}

Block* BlocksSystem::getBlockOfObject(Object* obj) {
	for(vector<Block*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		Object* temp = (*iter)->getElementByID(obj->getID());
		if(temp!=NULL)
			return (*iter);
	}
	return NULL;
}

void BlocksSystem::prepareBlocks() {
	for(vector<Block*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		Block* temp = (*iter);
		temp->calcCentreAndRadius();
	}
}

void BlocksSystem::buildMST() {

	set<Object*> allObjects;

	//collect the unique objects
	for(vector<Block*>::iterator iter = blocks.begin();iter!=blocks.end();iter++) {
		Block* temp = (*iter);
		set<Object*>* tempList = temp->getObjectsSet();
		allObjects.insert(tempList->begin(), tempList->end());
	}

	//create the disjoint set thing for the Kruskal algorithm
	/********************/
	/*std::vector<Element> elements_dsets_vect;
	elements_dsets_vect.reserve(allObjects.size());

	Rank rank(elements_dsets_vect);
	Parent parent(elements_dsets_vect);

	boost::disjoint_sets<Rank*, Parent*> sets(&rank, &parent);*/

	std::vector<int>  rank (Object::getMaxID()+1);
	std::vector<int>  parent (Object::getMaxID()+1);
	boost::disjoint_sets<int*,int*> ds(&rank[0], &parent[0]);

	/********************/

	GraphX* mst_shadow = new GraphX(1,0);
	GraphX* mst = new GraphX(1,0);
	for(set<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {
		Object* obj = (*iter);
		mst->AddNode(obj->getID(), 1, obj);
		mst_shadow->AddNode(obj->getID(), 1, obj);

		// for dsets
		ds.make_set(obj->getID());
	}

	vector<PossibleEdge*> finalAddedEdgesSorted;

	//create the minimum spanning tree
	int k = 1;
	while(true) {

		cout<<"k = "<<k<<endl;

		//sort objects based on the distance ascendingly
		vector<PossibleEdge*> allObjectsSorted;

		cout<<"Get distances ..."<<flush;

		int count = 0;

		for(set<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {
			Object* obj = (*iter);
			int coreK = 3;
			//pair<Object*, double>* p_o = obj->getKShortestRechability(NULL, this, coreK, k);
			vector<pair<Object*, double>*>* v_p_o = obj->getObjectsWithDistance(NULL, this, (k-1)/1.0, k/1.0);

			//put in place
			for(vector<pair<Object*, double>*>::iterator iter1 = v_p_o->begin();iter1!=v_p_o->end();iter1++) {

				pair<Object*, double>* p_o = (*iter1);

				vector<PossibleEdge*>::iterator innerIter = allObjectsSorted.begin();
				for(;innerIter!=allObjectsSorted.end();innerIter++) {
					if(p_o->second<(*innerIter)->dist)
						break;
				}

				PossibleEdge* pe = new PossibleEdge(obj, p_o->first, p_o->second);
				allObjectsSorted.insert(innerIter, pe);
			}
			count++;
//			cout<<count<<endl;
		}



		cout<<"Number of new edges = "<<": "<<allObjectsSorted.size()<<endl;
		cout<<"DONE."<<flush;

		cout<<"Generate MST ..."<<flush;
		//generate the mst tree
		vector<PossibleEdge*>::iterator innerIter = allObjectsSorted.begin();
		for(;innerIter!=allObjectsSorted.end();innerIter++) {
//			(*innerIter)->print(cout);

			Object* src = (*innerIter)->src;
			NodeX* srcNode = mst_shadow->getNodeWithID(src->getID());
			Object* dest = (*innerIter)->dest;
			NodeX* destNode = mst_shadow->getNodeWithID(dest->getID());
			double dist = (*innerIter)->dist;

			if(ds.find_set(src->getID())==ds.find_set(dest->getID()))
			{
				continue;
			}

			if(srcNode->isItConnectedWithNodeID(destNode->getID())) {
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

			mst_shadow->addEdge(srcNode->getID(), destNode->getID(), dist);
		}

		cout<<"DONE."<<flush;

		//check if it is connected
		if(mst->isConnected()) {
			cout<<*mst;
			break;
		}

//		cout<<"SHADOW SMT:"<<*mst_shadow<<endl;
//		cout<<"SMT:"<<*mst<<endl;

		k++;
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

	for(set<Object*>::iterator iter = allObjects.begin();iter!=allObjects.end();iter++) {
		Object* obj = (*iter);
		ds1.make_set(obj->getID());
	}

	for(vector<PossibleEdge*>::iterator iter1 = finalAddedEdgesSorted.begin();iter1!=finalAddedEdgesSorted.end();iter1++)
	{
		Object* src = (*iter1)->src;
		NodeX* srcNode = mst_shadow->getNodeWithID(src->getID());
		Object* dest = (*iter1)->dest;
		NodeX* destNode = mst_shadow->getNodeWithID(dest->getID());
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

		cout<<"*** "<<parent<<","<<hnSrc<<","<<hnDest<<endl;

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
	root->print();

	cout<<"Condensed advanced tree ..."<<endl;//<--- you have somthing wrong in your logic here!
	CH_Node* chn = new CH_Node();
	chn->setStart(0);
	chn->setWidth(root->getNumItems());
	root->condense(2, chn);
	chn->traverseAreas();
	chn->print(0);

	cout<<"Clusters ...."<<endl;
	vector<CH_Node*> clusters;
	chn->getClusters(&clusters);
	for(vector<CH_Node*>::iterator iter = clusters.begin();iter!=clusters.end();iter++) {
		(*iter)->print(0, true);
	}
}
