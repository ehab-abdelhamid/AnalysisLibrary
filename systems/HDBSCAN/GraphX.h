/*
 * GraphX.h
 *
 *  Created on: Mar 13, 2013
 *      Author: ehab
 *  This is for graph interface
 */

#ifndef GRAPHX_H_
#define GRAPHX_H_

#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <set>
#include <limits>
#include "NodeX.h"

#define allowThreaded 1

using namespace std;

class GraphX
{
private:
	int id;
	int type;//0-undirected, 1-directed, default value is 0. IMPORTANT: current algorithms do not support directed graphs
	map<int, NodeX*> nodes;
	//tr1::unordered_map<int, void* > nodesFeatures;
	int numOfEdges;
	tr1::unordered_map<double, set<int>* > nodeIDsByLabel;
	tr1::unordered_map<double, set<NodeX*>* > nodesByLabel;
	int maxLabel = -1;
	int maxNodeID = -1;


public:
	GraphX();
	GraphX(int, int );
	void init(int, int );
	GraphX(GraphX* );
	~GraphX();
	bool parseData(istream& data);
	bool loadFromFile(string fileName);
	bool loadFromString(string data);
	void loadNodesStats(string fileName);
	NodeX* AddNode(int id, double label, Object* );
	void addEdge(int , int , double);
	void addEdge(int , int , double ,map<string, void* >& );
	void removeEdge(int , int );
	void removeNode_IgnoreEdges(int );
	map<int,NodeX*>::const_iterator getNodesIterator();
	map<int,NodeX*>::const_iterator getNodesEndIterator();
	int getID() {return id;}
	int getType() {return type;}
	double getEdgeLabel(int srcNode, int destNode);
	int getNumOfNodes();
	NodeX* getNodeWithID(int nodeID);
	set<int>* getNodeIDsByLabel(double );
	set<NodeX*>* getNodesByLabel(double );
	friend ostream& operator<<(ostream& os, const GraphX& g);
	int getNumOfEdges() {return numOfEdges;}
	bool isConnected();

	int getNumDistinctLabels() { return nodeIDsByLabel.size(); }
};

ostream& operator<<(ostream& , const GraphX& );

#endif /* GRAPHX_H_ */
