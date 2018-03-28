/*
 * RandomBlock.cpp
 *
 *  Created on: 27 Nov 2017
 *      Author: eabdelha
 */

#include "RandomBlock.h"
#include "PossibleEdge.h"
#include "GraphX.h"
#include "TimeUtility.h"

struct less_than_key
{
    inline bool operator() (PossibleEdge*& struct1, PossibleEdge*& struct2)
    {
        return (struct1->dist < struct2->dist);
    }
};

struct less_than_key_object
{
    inline bool operator() (Object*& struct1, Object*& struct2)
    {
        return (struct1->getCoreDistance() < struct2->getCoreDistance());
    }
};

RandomBlock::RandomBlock() {
}

void RandomBlock::addElement(Object* obj) {
	elements.push_back(obj);
}

void RandomBlock::addObject(Object* obj) {
	objects.push_back(obj);
}

//void RandomBlock::partition(long r) {
//	partitions = ClosePartition::generatePartitions(objects, r);
//}

void RandomBlock::prepare() {

	//create the kdd tree
	int maxPts = objects.size();
	int dim = objects.at(0)->getValuesSize();

	dataPts1 = annAllocPts(maxPts, dim);
	int count = 0;
	for(vector<Object*>::iterator iter = objects.begin();iter!=objects.end();iter++) {
		ANNpoint tempPoint = dataPts1[count];
		for(int i=0;i<(*iter)->getValuesSize();i++) {
			tempPoint[i] = (*iter)->getValueByIndex(i);
		}
		count++;
	}
	kdTreeObjects = new ANNkd_tree(dataPts1, maxPts, dim);

	//create the kdd tree
	maxPts = elements.size();
	dim = elements.at(0)->getValuesSize();

	dataPts2 = annAllocPts(elements.size(), dim);
	count = 0;
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		ANNpoint tempPoint = dataPts2[count];
		for(int i=0;i<(*iter)->getValuesSize();i++) {
			tempPoint[i] = (*iter)->getValueByIndex(i);
		}
		count++;
	}
	kdTreeElements = new ANNkd_tree(dataPts2, maxPts, dim);
}

void RandomBlock::prepareSortedElements() {

	//create the sorted list
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		sortedElements.push_back((*iter));
	}

	//sort
	std::sort(sortedElements.begin(), sortedElements.end(), less_than_key_object());
}

void RandomBlock::print(ostream& out) {
//	out<<"Core elements = "<<this->elements.size()<<", Others = "<<this->objects.size()<<", #partitions = "<<partitions.size()<<endl;
	out<<"Core elements = "<<this->elements.size()<<", Others = "<<this->objects.size()<<endl;
//	for(vector<ClosePartition*>::iterator iter = partitions.begin();iter!=partitions.end();iter++) {
//		cout<<"\t";
//		(*iter)->print(cout);
//	}
}

/*vector<PossibleEdge*> RandomBlock::getPairsWithDistBound(double min, double max) {
	vector<PossibleEdge*> allElements;

	int total = 0;
	int pruned = 0;
	int totalEdges = 0;
	int remainingEdges = 0;
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		for(vector<ClosePartition*>::iterator iter1 = partitions.begin();iter1!=partitions.end();iter1++) {
			ClosePartition* cp = (*iter1);
			double dist = cp->getDistanceToCenter(*iter);
			if(max<dist-cp->getR() || min>dist+cp->getR()) {
				//the case when the iter element cannot be within a distance with any element in cp
				pruned++;
			} else {
				//now search for candidate pairs
				vector<PossibleEdge*> partialList = cp->getElementsWithDistBound(*iter, min, max);

				totalEdges+=partialList.size();
				//remove edges that connect objects within the same connected component
				for(vector<PossibleEdge*>::iterator pruningIter = partialList.begin();pruningIter!=partialList.end();) {
					PossibleEdge* pe = (*pruningIter);
					if(ds->find_set(pe->src->getID())==ds->find_set(pe->dest->getID()))
					{
						pruningIter = partialList.erase(pruningIter);
					}
					else
					{
						ds->union_set(pe->src->getID(), pe->dest->getID());
						pruningIter++;
					}
				}
				remainingEdges+=partialList.size();

				//apply merge sort
				vector<PossibleEdge*>::iterator allElementsIter = allElements.begin();
				vector<PossibleEdge*>::iterator partialListIter = partialList.begin();
				while(true) {
					if(allElementsIter==allElements.end()) {
						allElements.insert(allElementsIter, partialListIter, partialList.end());
						break;
					}
					if(partialListIter==partialList.end()) {
						break;
					}

					if((*allElementsIter)->dist<(*partialListIter)->dist) {
						allElementsIter++;
					} else {
						allElementsIter = allElements.insert(allElementsIter, (*partialListIter));
						allElementsIter++;
						partialListIter++;
					}
				}

				//std::copy(partialList.begin(), partialList.end(), std::back_insert_iterator<std::vector<PossibleEdge*> >(allElements));
			}
			total++;
		}
	}

	cout<<"(pruned/total)="<<pruned<<"<"<<total<<"("<<((double(pruned))/total)<<")"<<endl;
	cout<<"(remainingEdges/totalEdges)="<<remainingEdges<<"<"<<totalEdges<<"("<<((double(remainingEdges))/totalEdges)<<")"<<endl;

	return allElements;
}*/

void RandomBlock::computeCoreDistances(int nCorePoints) {
	ANNidxArray nnIdx;
	ANNdistArray dists;
	int k = nCorePoints+1;//considering the point itself, we return nCorePoints+1
	nnIdx = new ANNidx[k];
	dists = new ANNdist[k];

	for(vector<Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		Object* object = *it;

		kdTreeObjects->annkSearch(object->getValuesAsANNPoint(), k, nnIdx, dists, 0);

		if(dists[0]!=0)
			cout<<"ERRRORRRRRR284728432: the first point, I thought, is supposed to be the point itself!"<<endl;

		double maxDistance = dists[k-1];
		/*for(int i=0;i<k;i++) {
			if(dists[i]>maxDistance)
				maxDistance = dists[i];
			else
				cout<<"ERR98921121121"<<endl;
		}*/

		object->setCoreDistance(sqrt(maxDistance));
	}

	delete[] nnIdx;
	delete[] dists;
}

void RandomBlock::getCoreDistance(Object* obj, int knn, vector<pair<Object*,double>* >* sortedList) {

	for(std::vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it) {
		Object* otherObject = *it;
		if(obj==otherObject)
			continue;

		double distance = obj->distanceTo(otherObject);
		//put it in order
		int count = 0;
		bool insert = true;
		std::vector<pair<Object*,double>*>::iterator its = sortedList->begin();
		for(; its != sortedList->end(); ++its) {
			pair<Object*, double>* o = *its;
			if(o->first==otherObject) {
				insert = false;
				break;
			}
			if(o->second>distance)
				break;
			count++;
			if(knn!=-1 && count>knn)
				break;
		}

		if(insert)
			sortedList->insert(its, new pair<Object*, double>(otherObject, distance));

		if(knn!=-1 && sortedList->size()>knn)
		{
			delete sortedList->back();
			sortedList->pop_back();
		}
	}
}

void RandomBlock::buildLocalMST() {

	vector<PossibleEdge*> orderedEdges;

	for(std::vector<Object*>::iterator it1 = elements.begin(); it1 != elements.end(); ++it1) {
		for(std::vector<Object*>::iterator it2 = elements.begin(); it2 != elements.end(); ++it2) {
			if((*it1)==(*it2))
				continue;

			PossibleEdge* pe = new PossibleEdge((*it1), (*it2), (*it1)->distanceTo((*it2)));
			orderedEdges.push_back(pe);
		}
	}

	//sort them
	std::sort(orderedEdges.begin(), orderedEdges.end(), less_than_key());
//	for(vector<PossibleEdge*>::iterator iter = orderedEdges.begin();iter!=orderedEdges.end();iter++) {
//		(*iter)->print(cout);
//	}

	//build MST
	std::vector<int>  rank (objects.size()+1);
	std::vector<int>  parent (objects.size()+1);
	boost::disjoint_sets<int*,int*> ds(&rank[0], &parent[0]);

	////////////////////////////////////

	GraphX* mst = new GraphX(1,0);
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		Object* obj = (*iter);
		mst->AddNode(obj->getID(), 1, obj);

		// for dsets
		ds.make_set(obj->getID());
	}

	while(true) {
		//pick the edge with smallest distance
		PossibleEdge* pe = *(orderedEdges.begin());
		orderedEdges.erase(orderedEdges.begin());

		//try to add this edge if possible
		if(ds.find_set(pe->src->getID())!=ds.find_set(pe->dest->getID())) {
			mst->addEdge(pe->src->getID(), pe->dest->getID(), pe->dist);
			ds.union_set(pe->src->getID(), pe->dest->getID());
		}

		if(orderedEdges.size()%1000==0) {

			for(vector<PossibleEdge*>::iterator iter = orderedEdges.begin();iter!=orderedEdges.end();) {
				if(ds.find_set((*iter)->src->getID())==ds.find_set((*iter)->dest->getID()))
					iter = orderedEdges.erase(iter);
				else
					iter++;
			}

//			cout<<orderedEdges.size()<<", "<<mst->getNumOfNodes()<<", "<<mst->getNumOfEdges()<<endl;
		}

		//break if all nodes are included in mst
		if(mst->getNumOfEdges()>=(mst->getNumOfNodes()-1) && mst->isConnected())
			break;
	}
}

vector<PossibleEdge*>* RandomBlock::buildLocalMST_Prims() {

	//generate a list of all to be checked nodes
	vector<PossibleEdge*> candidates;//implement candidates as a set tehn apply sorting on this list
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		candidates.push_back(new PossibleEdge((*iter), NULL, -1));
	}

//	GraphX mst = new GraphX();

	//start with the first element
	PossibleEdge* currElement = *(candidates.begin());
	candidates.erase(candidates.begin());
//	mst.AddNode(currElement->src->getID(), 1);

	vector<PossibleEdge*>* orderedEdges = new vector<PossibleEdge*>();

	while(true) {

		//make sure the candidates list is always sorted
		/*vector<PossibleEdge*>::iterator iter2 = candidates.begin();
		double lastDist = (*iter2)->dist;
		iter2++;
		for(;iter2!=candidates.end();iter2++) {
			//cout<<"U: "<<lastDist<<", "<<(*iter2)->dist<<endl<<flush;
			if(lastDist>(*iter2)->dist) {
				cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB: "<<lastDist<<", "<<(*iter2)->dist<<endl<<flush;
				return NULL;
			}
			lastDist = (*iter2)->dist;
		}*/

		//calculate distances from the newly added node to the list of candidates
//		for(vector<PossibleEdge*>::iterator iter = candidates.begin();iter!=candidates.end();iter++) {

		//for(int i=candidates.size()-1;i>=0;i--) {
		for(int i=0;i<candidates.size();i++) {
			PossibleEdge* pe = candidates[i];

			//there is no need to check distance if pe distance is less than the core distance of currElement (i.e., no way to decrease pe->distance)
			if(pe->dist!=-1 && pe->dist<currElement->src->getCoreDistance())
				continue;

			double dist = currElement->src->coreDistanceTo(pe->src, pe->dist);

			if(pe->dist==-1 || pe->dist>dist) {

				//update the distance of current element
				pe->dest = currElement->src;
				pe->dist = dist;

				if(i==0) continue;

				//update its position
				int j = 0;

				int start = 0;
				int end = i;
				int mid = 0;
				while(true) {
					mid = start + (end-start+1)/2;
					if(end-start<2) {
						if(candidates[start]->dist>dist)
							j = start;
						else
							j = end;
						break;
					}
					if(dist>candidates[mid]->dist) {
						start = mid;
					} else {
						end = mid;
					}
				}

//				for(int k=i-1;k>=j;k--)
//					candidates[k+1] = candidates[k];
				std::memmove(candidates.data()+(j+1), candidates.data()+j, sizeof(PossibleEdge*)*(i-j));
				candidates[j] = pe;

//				for(int l=1;l<candidates.size();l++) {
//					if(candidates[l]->dist==-1)
//						break;
//					if(candidates[l-1]->dist>candidates[l]->dist)
//					{
//						cout<<"FFFFFFFFFF "<<dist<<", l = "<<l<<", "<<candidates[l-1]->dist<<","<<candidates[l]->dist<<endl;
//						return NULL;
//					}
//				}
			}
		}

		delete currElement;
		//pick the first element in the list

		currElement = *(candidates.begin());
		candidates.erase(candidates.begin());
//		mst.AddNode(currElement->src->getID(), 1);
//		mst.addEdge(currElement->src->getID(), currElement->dest->getID(), currElement->dist);

		//add the edge ordered
		PossibleEdge* nePE = new PossibleEdge(currElement->src, currElement->dest, currElement->dist);
		/*list<PossibleEdge*>::iterator orderingIter = orderedEdges->begin();
		for(;orderingIter!=orderedEdges->end();orderingIter++) {
			if(nePE->dist<(*orderingIter)->dist) {
				break;
			}
		}
		orderedEdges->insert(orderingIter, nePE);
		*/
		orderedEdges->push_back(nePE);

		if(candidates.size()<=0) {
			delete currElement;
			break;
		}
	}

	std::sort(orderedEdges->begin(), orderedEdges->end(), less_than_key());

	return orderedEdges;
}

/**
 * This function is optimized by:
 * 1- KD tree for finding the closest object
 * 2- optimization#1: if objCoreDistance>=distance && objCoreDistance>=otherObjectCoreDistance then distance = objCoreDistance (i.e., objCoreDistance is the minimum possible distance)
 */
//with ANN tree index
std::pair<Object*, double>* RandomBlock::getClosestElement_CoreDistance(Object* obj) {
	int k = 1;
//	ANNidxArray nnIdx;
//	ANNdistArray dists;
//	nnIdx = new ANNidx[k];
//	dists = new ANNdist[k];
	ANNidx nnIdx;
	ANNdist dists;

	kdTreeElements->annkSearch(obj->getValuesAsANNPoint(), k, &nnIdx, &dists, 0);

	double distance = sqrt(dists);//sqrt(dists[0]);
	Object* otherObject = elements[nnIdx];//elements[nnIdx[0]];

	double objCoreDistance = obj->getCoreDistance();
	double otherObjectCoreDistance = otherObject->getCoreDistance();

	if(objCoreDistance>=distance && objCoreDistance>=otherObjectCoreDistance) {

//		double min1 = std::numeric_limits<double>::max();
//		Object* otherObject1 = NULL;
//
//		for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
//
//			double tempDistance = obj->coreDistanceTo((*iter));
//			if(min1>tempDistance) {
//				min1 = tempDistance;
//				otherObject1 = (*iter);
//			}
//		}
//
//		if(objCoreDistance!=min1 || otherObject!=otherObject1)
//			cout<<"Optimized decision: "<<otherObject<<", "<<objCoreDistance<<"\n Slow Decision: "<<otherObject1<<", "<<min1<<endl;

		return new std::pair<Object*, double>(otherObject, objCoreDistance);
	}
	else {
		double min = distance;
		if(min<otherObjectCoreDistance) min = otherObjectCoreDistance;
		Object* otherObject1 = otherObject;

		for(vector<Object*>::iterator iter = sortedElements.begin();iter!=sortedElements.end();iter++) {

			if(min<(*iter)->getCoreDistance())
				break;

			if(objCoreDistance>(*iter)->getCoreDistance())
				continue;

			double tempDistance = obj->coreDistanceTo((*iter), min);
			if(min>tempDistance) {
				min = tempDistance;
				otherObject1 = (*iter);
			}
		}

		return new std::pair<Object*, double>(otherObject1, min);
	}
}
//with no ANN tree index
/*std::pair<Object*, double>* RandomBlock::getClosestElement_CoreDistance(Object* obj) {
	int k = 1;

	double distance = 100000000;
	Object* otherObject = elements[0];

	double objCoreDistance = obj->getCoreDistance();
	double otherObjectCoreDistance = otherObject->getCoreDistance();

	double min = distance;
	if(min<otherObjectCoreDistance) min = otherObjectCoreDistance;
	Object* otherObject1 = otherObject;

	for(vector<Object*>::iterator iter = sortedElements.begin();iter!=sortedElements.end();iter++) {

		if(min<(*iter)->getCoreDistance())
			break;

		if(objCoreDistance>(*iter)->getCoreDistance())
			continue;

		double tempDistance = obj->coreDistanceTo((*iter), min);
		if(min>tempDistance) {
			min = tempDistance;
			otherObject1 = (*iter);
		}
	}

	return new std::pair<Object*, double>(otherObject1, min);
}*/

RandomBlock::~RandomBlock() {
	delete kdTreeObjects;
	delete kdTreeElements;

	//pa[i] = &(p[i*dim]);
	delete[] *dataPts1;
	delete[] *dataPts2;

	delete[] dataPts1;
	delete[] dataPts2;
}
