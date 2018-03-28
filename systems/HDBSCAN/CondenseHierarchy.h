/*
 * CondenseHierarchy.h
 *
 *  Created on: 30 Oct 2017
 *      Author: eabdelha
 */

#ifndef CONDENSEHIERARCHY_H_
#define CONDENSEHIERARCHY_H_

#include<vector>
#include "NodeX.h"

using namespace std;

class H_Node;

class CH_Node {
private:
	double start;
	double end;
	int width;
	CH_Node* component = NULL;
	vector<CH_Node*> children;
	H_Node* dataNode;
	double tempArea = -1;

public:

	CH_Node() { dataNode = NULL; }
	~CH_Node() { }
	void setStart(double start);
	void setEnd(double end);
	void setWidth(int width);
	void setDataNode(H_Node* hnode);
	H_Node* getDataNode();
	int getWidth();
	double getArea();
	void addChild(CH_Node* node);
	void setComponent(CH_Node* node);
	void print(int depth, bool single=false, bool isComponent=false);
	void printAsList(bool single=false, bool isComponent=false);
	double traverseAreas();
	void getClusters(vector<CH_Node*>* clusters);
	static void deleteCh_Node(CH_Node* node);
};


#endif /* CONDENSEHIERARCHY_H_ */
