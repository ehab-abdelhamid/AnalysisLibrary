/*
 * ClosePartition.cpp
 *
 *  Created on: 27 Nov 2017
 *      Author: eabdelha
 */

#include "ClosePartition.h"
#include"PossibleEdge.h"

vector<ClosePartition*> ClosePartition::generatePartitions(vector<Object*> outObjs, double r) {
	vector<ClosePartition*> partitions;

	int count = 0;
	for(vector<Object*>::iterator iter = outObjs.begin();iter!=outObjs.end();iter++) {
		bool added = false;
//		cout<<count<<endl;
		count++;
		Object* obj = (*iter);
		for(vector<ClosePartition*>::iterator iter1 = partitions.begin();iter1!=partitions.end();iter1++) {
			if((*iter1)->addObject(obj)) {
				added = true;
				break;
			}
		}
		if(!added) {
			ClosePartition* cp = new ClosePartition();
			cp->setR(r);
			cp->addObject(obj);
			partitions.push_back(cp);
		}
	}

	//make sure your partitioning technique is correct (just for debug)
	int sum = 0;
	for(vector<ClosePartition*>::iterator iter = partitions.begin();iter!=partitions.end();iter++) {
		sum+=(*iter)->objects.size();
	}
//	cout<<"#allObjects = "<<outObjs.size()<<"#partitions = "<<partitions.size()<<", sum = "<<sum<<endl;

	return partitions;
}

vector<PossibleEdge*> ClosePartition::getElementsWithDistBound(Object* obj, double min, double max) {
	vector<PossibleEdge*> r;

	for(vector<Object*>::iterator iter = objects.begin();iter!=objects.end();iter++) {
		long dist = obj->distanceTo((*iter));
		if(dist>=min && dist<max) {
			PossibleEdge* pe = new PossibleEdge(obj, (*iter), dist);

			vector<PossibleEdge*>::iterator sortingIter = r.begin();
			for(;sortingIter!=r.end();sortingIter++) {
				if(pe->dist>(*sortingIter)->dist)
					break;
			}
			r.insert(sortingIter, pe);
//			r.push_back(pe);
		}
	}

	return r;
}

double ClosePartition::getDistanceToCenter(Object* obj) {
	return centroid->distanceTo(obj);
}

bool ClosePartition::isWithinDist(Object* obj) {
	return centroid->distanceTo(obj)<r;
}

bool ClosePartition::addObject(Object* obj) {
	if(objects.size()==0) {
		centroid = obj;
		objects.push_back(obj);
		return true;
	}

	if(!isWithinDist(obj))
		return false;

	objects.push_back(obj);
	return true;
}

void ClosePartition::print(ostream& out) {
	centroid->print(out);
	out<<", "<<objects.size()<<", "<<r;
	out<<endl;
}


