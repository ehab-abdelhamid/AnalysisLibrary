/*
 * Hierarchy.h
 *
 *  Created on: 24 Oct 2017
 *      Author: eabdelha
 */

#ifndef HIERARCHY_H_
#define HIERARCHY_H_

#include<vector>
#include<set>
#include "NodeX.h"
#include "CondenseHierarchy.h"

using namespace std;

class H_Node;

class HNode_Distance {
public:
	HNode_Distance(H_Node* hnode, double distance) {this->hnode = hnode; this->distance = distance;}
	H_Node* hnode;
	double distance;
};

class H_Node {
private:
	vector<HNode_Distance*> elements;
	NodeX* nodeValue = NULL;
	int numItems = -1;//number items (not elements) included in all of its sub-clusters
	//set<NodeX*> items;
	H_Node* parent = NULL;//initially parent is not set, it is set in the setNumItems method

public:
	~H_Node();
	int getNumItems();
	set<NodeX*>* getItems();
	void deleteDescendents();
	void setNumItems(H_Node* parent = NULL);
	void addNode(H_Node* hn, double distance);
	void setNode(NodeX* n);
	int getNodeID();
	void condense(int minSize, CH_Node* ch, H_Node* parent = NULL);
	void print(H_Node* parent = NULL);
	void printAsList(H_Node* parent = NULL);
};



#endif /* HIERARCHY_H_ */
