/*
 * GraphX.cpp
 *
 *  Created on: Mar 13, 2013
 *      Author: ehab
 */

#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include<iostream>
#include<fstream>
#include<sstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include "GraphX.h"
#include "EdgeX.h"
#include "Utils.h"

/**
 * Constructor
 */
GraphX::GraphX()
{
	init(0,0);
}

GraphX::GraphX(int id, int type)
{
	init(id, type);
}

void GraphX::init(int id, int type)
{
	this->id = id;
	this->type = type;
	if(type==1)
	{
		cout<<"IMPORTANT! current set of algorithms do not support directed graphs!";
		exit(0);
	}
	numOfEdges = 0;
}

GraphX::GraphX(GraphX* graph)
{
	this->id = graph->getID();
	this->type = graph->getType();
	numOfEdges = graph->getNumOfEdges();

	//copy nodes
	for(map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		NodeX* oldNode = iter->second;
		this->AddNode(oldNode->getID(), oldNode->getLabel(), oldNode->getObject());
	}

	//copy edges
	for(map<int, NodeX*>::const_iterator iter = graph->getNodesIterator();iter!=graph->getNodesEndIterator();++iter)
	{
		NodeX* oldNode = iter->second;
		NodeX* node = nodes.at(oldNode->getID());

		for(map<int, void*>::iterator iter2 = oldNode->getEdgesIterator();iter2!=oldNode->getEdgesEndIterator();++iter2)
		{
			EdgeX* edge = (EdgeX*)(iter2->second);
			node->addEdge(nodes.at(edge->getOtherNode()->getID()) , edge->getLabel(), this->type);
		}
	}
}

/**
 * Add a node to the graph
 * Parameters are: node id, and node label
 */
NodeX* GraphX::AddNode(int id, double label, Object* obj)
{
	//assign the maximum label value (note: the current node label implementation is double!)
	if(maxLabel<label)
		maxLabel = label;

	if(maxNodeID<id)
		maxNodeID = id;

	map<int, NodeX*>::iterator temp = nodes.find(id);
	if(temp!=nodes.end())
		return temp->second;

	NodeX* node = new NodeX(id, label);
	node->setObject(obj);
	nodes.insert(std::pair<int, NodeX*>(id, node));

	//add to the 'node IDs by label' map
	tr1::unordered_map<double, set<int>* >::iterator iter = nodeIDsByLabel.find(label);
	if(iter==nodeIDsByLabel.end())
	{
		nodeIDsByLabel.insert(std::pair<double, set<int>*>(label, new set<int>()));
		iter = nodeIDsByLabel.find(label);
	}
	iter->second->insert(id);

	//add to the 'nodes by label' map
	tr1::unordered_map<double, set<NodeX*>* >::iterator iter1 = nodesByLabel.find(label);
	if(iter1==nodesByLabel.end())
	{
		nodesByLabel.insert(std::pair<double, set<NodeX*>*>(label, new set<NodeX*>()));
		iter1 = nodesByLabel.find(label);
	}
	iter1->second->insert(node);

	return node;
}

/**
 * Add an edge between two nodes in the graph, the nodes must exist before adding the edge
 * Parameters: the source node ID, the destination node ID, and the edge label
 * For undirected graphs, one more edge will be added in the reverse direction
 */
void GraphX::addEdge(int srcID, int destID, double edgeLabel)
{
//	cout<<"src = "<<srcID<<"' dest = "<<destID<<", dist = "<<edgeLabel<<endl<<flush;
	nodes[srcID]->addEdge(nodes[destID], edgeLabel, this->type);

	if(this->type==0)
	{
		nodes[destID]->addEdge(nodes[srcID], edgeLabel, this->type);
	}

	numOfEdges++;
}

void GraphX::removeEdge(int id1, int id2)
{
	nodes[id1]->removeEdge(nodes[id2], this->type);
	if(nodes[id1]->getEdgesSize()==0)
		removeNode_IgnoreEdges(id1);

	if(this->type==0)
	{
		nodes[id2]->removeEdge(nodes[id1], this->type);
		if(nodes[id2]->getEdgesSize()==0)
			removeNode_IgnoreEdges(id2);
	}

	numOfEdges--;
}

/**
 * remove a node from the graph, ignoring edges (as a prepost all edges connecting to this node should be already removed)
 */
void GraphX::removeNode_IgnoreEdges(int nodeID)
{
	nodeIDsByLabel.find(getNodeWithID(nodeID)->getLabel())->second->erase(nodeID);
	delete nodes.find(nodeID)->second;
	nodes.erase(nodeID);
}

void GraphX::addEdge(int srcID, int destID, double edgeLabel, map<string, void* >& edgeToFreq)
{
	//cout<<"Adding edge: "<<srcID<<"---"<<destID<<"\n";
	this->addEdge(srcID, destID, edgeLabel);
}

/**
 * Load a graph file that has .lg format
 * return true if loading is done correctly, otherwise false
 */
bool GraphX::loadFromFile(string fileName)
{
	cout<<"Loading graph from file: "<<fileName<<endl;
	ifstream file (fileName.c_str(), ios::in);
	if(!file)
	{
		cout << "While opening a file an error is encountered" << endl;
		return false;
    }

	if(!parseData(file))
		return false;

	file.close();

	return true;
}

bool GraphX::loadFromString(string data)
{
	cout<<"Loading graph from string: "<<data<<endl;
	istringstream str(data);

	bool b = parseData(str);

	if(!b)
		return false;

	return true;
}

bool GraphX::parseData(istream& data)
{
	//read the first line
	char temp_ch;
	data>>temp_ch;data>>temp_ch;data>>temp_ch;

	int numEdgesLoaded = 0;

	while (true)
	{
		char ch;
		data>>ch;

		if(data.eof()) break;

		//to add nodes
		if(ch=='v')
		{
			int id;
			double label;
			data>>id;
			data>>label;
//			if(id<100000)//remove me
				this->AddNode(id, label, NULL);
		}
		else if(ch=='e')//to add edges
		{
			int id1;
			int id2;
			double label;
			data>>id1;
			data>>id2;
			data>>label;
//			if(id1<100000 && id2<100000)//remove me
			{
			this->addEdge(id1, id2, label);
			if(numEdgesLoaded%1000000==0)
				cout<<"Loaded edges = "<<numEdgesLoaded<<endl;
			numEdgesLoaded++;
			}
		}
	}

	return true;
}

NodeX* GraphX::getNodeWithID(int nodeID)
{
	map<int, NodeX*>::iterator iter = nodes.find(nodeID);
	if(iter==nodes.end())
		return NULL;
	else
		return iter->second;
}

set<int>* GraphX::getNodeIDsByLabel(double label)
{
	tr1::unordered_map<double, set<int>* >::iterator iter = nodeIDsByLabel.find(label);
	if(iter==nodeIDsByLabel.end())
		return NULL;
	return iter->second;
}

set<NodeX*>* GraphX::getNodesByLabel(double label)
{
	tr1::unordered_map<double, set<NodeX*>* >::iterator iter = nodesByLabel.find(label);
	if(iter==nodesByLabel.end())
		return NULL;
	return iter->second;
}

map<int,NodeX*>::const_iterator GraphX::getNodesIterator()
{
	return nodes.begin();
}

map<int,NodeX*>::const_iterator GraphX::getNodesEndIterator()
{
	return nodes.end();
}

/**
 * get edge label from the scr node to the destination node
 * if either the src node is not found, or the detination node is not found from the src node, then return 0.0001
 */
double GraphX::getEdgeLabel(int srcNodeID, int destNodeID)
{
	NodeX* srcNode;
	map<int, NodeX* >::iterator temp = nodes.find(srcNodeID);
	if(temp==nodes.end())
			return 0.0001;
	srcNode = (*temp).second;

	EdgeX* edge = (EdgeX*)srcNode->getEdgeForDestNode(destNodeID);
	if(edge==NULL)
		return 0.0001;

	return edge->getLabel();
}

int GraphX::getNumOfNodes()
{
	return this->nodes.size();
}

bool GraphX::isConnected()
{
	if(nodes.size()<2)
		return true;

	//start from any node
	NodeX* node = nodes.begin()->second;
	map<int,NodeX*> visited;
	map<int,NodeX*> toVisit;

	toVisit.insert(std::pair<int, NodeX*>(node->getID(), node));
	while(toVisit.size()>0)
	{
		//1- pop a node for the to be visited list, 2- remove it from toVisit, and 3- add it to the visited list
		node = toVisit.begin()->second;//1
		toVisit.erase(toVisit.begin());//2
		visited.insert(std::pair<int, NodeX*>(node->getID(), node));//3
		//add its neighbours
		for(map<int, void*>::iterator iter = node->getEdgesIterator();iter!=node->getEdgesEndIterator();++iter)
		{
			int id = iter->first;
			if(visited.find(id)!=visited.end())
				continue;
			EdgeX* edge = (EdgeX*)iter->second;
			toVisit.insert(std::pair<int, NodeX*>(id, edge->getOtherNode()));
		}
	}

	if(visited.size()==nodes.size())
		return true;
	else {
//		cout<<"Graph is disconnected. The connected part covers ";
		//cout<<(((double)visited.size())/nodes.size())<<endl;
//		cout<<(((double)this->numOfEdges)/(nodes.size()-1))<<endl;
		return false;
	}

}

/**
 * get the canonical label of this graph
 */
//string GraphX::getCanonicalLabel()
//{
//	if(CL.length()>0)
//		return CL;
//	CL = CanonicalLabel::generate(this);
//	return CL;
//}

//output graph
//ostream& operator<<(ostream& os, const GraphX& g)
//{
//	for( map<int,NodeX*>::const_iterator ii=g.nodes.begin(); ii!=g.nodes.end(); ++ii)
//	{
//		os<<*((*ii).second)<<endl;
//	}
//    return os;
//}

ostream& operator<<(ostream& os, const GraphX& g)
{
	os<<"t # 1\n";
	//output the nodes
	for( map<int,NodeX*>::const_iterator ii=g.nodes.begin(); ii!=g.nodes.end(); ++ii)
	{
		NodeX* node = ii->second;
		os<<"v "<<node->getID()<<" "<<node->getLabel()<<"\n";
	}

	//output the edges
	tr1::unordered_set<string> savedEdges;//list to keep track of already saved edges
	for( map<int,NodeX*>::const_iterator ii=g.nodes.begin(); ii!=g.nodes.end(); ++ii)
	{
		NodeX* node = ii->second;
		for(map<int, void*>::iterator iter1 = node->getEdgesIterator();iter1!=node->getEdgesEndIterator();iter1++)
		{
			EdgeX* edge = (EdgeX*)iter1->second;

			//check whether it has been added before or not
			string sig = intToString(node->getID())+"_"+doubleToString(edge->getLabel())+"_"+intToString(edge->getOtherNode()->getID());
			if(node->getID()<edge->getOtherNode()->getID())
				sig = intToString(edge->getOtherNode()->getID())+"_"+doubleToString(edge->getLabel())+"_"+intToString(node->getID());
			if(savedEdges.find(sig)==savedEdges.end())
			{
				savedEdges.insert(sig);
				os<<"e "<<node->getID()<<" "<<edge->getOtherNode()->getID()<<" "<<edge->getLabel()<<"\n";
			}
		}
	}

    return os;
}

//Destructor
GraphX::~GraphX()
{
	for( map<int,NodeX*>::iterator ii=nodes.begin(); ii!=nodes.end(); ++ii)
	{
		delete ii->second;
	}
	nodes.clear();

	for( tr1::unordered_map<double, set<int>* >::iterator ii=nodeIDsByLabel.begin(); ii!=nodeIDsByLabel.end(); ++ii)
	{
		delete ii->second;
	}
	nodeIDsByLabel.clear();

	for( tr1::unordered_map<double, set<NodeX*>* >::iterator ii=nodesByLabel.begin(); ii!=nodesByLabel.end(); ++ii)
	{
		delete ii->second;
	}
	nodesByLabel.clear();
}
