/*
 * NodeX.h
 *
 *  Created on: Mar 13, 2013
 *      Author: ehab
 *  This is for node interface
 */

#ifndef NODEX_H_
#define NODEX_H_

#include <map>
#include<vector>
#include<ostream>
#include "Object.h"

using namespace std;

class NodeX
{
private:
	int id;
	double label;
	string labelStr;
	map<int, void*> edges;
	map<int, void*> revEdges;
	Object* object;

public:
	NodeX(int id, double value);
	~NodeX();
	void addEdge(NodeX* , double edgeLabel, int graphType);
	void removeEdge(NodeX* , int graphType);
	friend ostream& operator<<(ostream& os, const NodeX& n);
	int getID() {return id;}
	double getLabel() {return label;}
	string getLabelStr() { return labelStr; }
	map<int, void*>::iterator getEdgesIterator() {return edges.begin();}
	map<int, void*>::iterator getEdgesEndIterator() {return edges.end();}
	int getEdgesSize() {return edges.size(); }
	void* getEdgeForDestNode(int destNodeID );
	bool isItConnectedWithNodeID(int nodeID);
	bool isItConnectedWithNodeID(int nodeID, double label);
	bool isNeighborhoodConsistent(NodeX* );
	void setObject(Object* obj) { this->object = obj; }
	Object* getObject() { return object; }
};

#endif /* NODEX_H_ */
