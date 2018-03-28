/*
 * Cluster.h
 *
 *  Created on: 9 Jan 2018
 *      Author: eabdelha
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include<vector>
#include "Hierarchy.h"

using namespace std;

class Cluster {
private:
	int cluserID = -1;
	vector<Object* > elements;

public:
	void setClusterID(int cID) { this->cluserID = cID; }
	int getClusterID() { return this->cluserID; }
	void addElement(Object* obj) {
		elements.push_back(obj);
	}

	void printElemValues(ostream& out) {
		for(vector<Object* >::iterator iter = elements.begin();iter!=elements.end();iter++) {
			out<<(*iter)->getID()<<",";
			(*iter)->printValues(out);
			out<<",";
		}
	}

	void populateObjects(H_Node* hnode) {
		set<NodeX*>* items = hnode->getItems();
		for(set<NodeX*>::iterator iter=items->begin();iter!=items->end();iter++) {
			NodeX* node = (*iter);
			addElement(node->getObject());
		}
		delete items;
	}

	bool exists(Object* obj) {
		//you can use a better implementation
		for(vector<Object* >::iterator iter = elements.begin();iter!=elements.end();iter++) {
			if(obj->getID()==(*iter)->getID())
				return true;
		}
		return false;
	}

	vector<Object* >::iterator getElementsBegin() {
		return elements.begin();
	}

	vector<Object* >::iterator getElementsEnd() {
		return elements.end();
	}

};


#endif /* CLUSTER_H_ */
