/*
 * CondenseHierarchy.cpp
 *
 *  Created on: 9 Jan 2018
 *      Author: eabdelha
 */

#include<iostream>
#include "CondenseHierarchy.h"
#include "Hierarchy.h"

using namespace std;

void CH_Node::setStart(double start) {
	this->start = start;
}

void CH_Node::setEnd(double end) {
	this->end = end;
}

void CH_Node::setWidth(int width) {
	this->width = width;
}

void CH_Node::setDataNode(H_Node* hnode) {
	dataNode = hnode;
}

H_Node* CH_Node::getDataNode() {
	return dataNode;
}

int CH_Node::getWidth() {
	return width;
}

double CH_Node::getArea() {
	double componentArea = 0;
	if(component!=NULL)
		componentArea = component->getArea();

	return (end-start)*width + componentArea;
}

void CH_Node::addChild(CH_Node* node) {
	children.push_back(node);
}

void CH_Node::setComponent(CH_Node* node) {
	this->component = node;
}

void CH_Node::print(int depth, bool single, bool isComponent) {
	for(int i=0;i<depth;i++)
		cout<<"\t";
	if(isComponent)
		cout<<"*";//to indicate this is a component
	cout<<"start = "<<start<<", end = "<<end<<", width = "<<width<<", area = "<<getArea()<<", tempArea = "<<tempArea<<endl;
	if(single) return;

	if(component!=NULL) {
		component->print(depth, false, true);
	}
	else {
		for(vector<CH_Node*>::iterator iter = children.begin();iter!=children.end();iter++) {
			(*iter)->print(depth+1);
		}
	}
}

void CH_Node::printAsList(bool single, bool isComponent) {

	if(dataNode!=NULL)
		dataNode->printAsList();
	cout<<endl;
	/*for(vector<NodeX*>::iterator iter = dataNodes.begin();iter!=dataNodes.end();iter++) {
		cout<<(*iter)->getID()<<", ";
	}

	if(single) return;

	if(component!=NULL) {
		component->printAsList(false, true);
	}
	else {
		for(vector<CH_Node*>::iterator iter = children.begin();iter!=children.end();iter++) {
			(*iter)->printAsList();
		}
	}*/
}

double CH_Node::traverseAreas() {
	double area = this->getArea();

	//get areas of the children
	double childrenArea = 0;
	CH_Node* currentCHNode = this;
	while(currentCHNode->component!=NULL)
		currentCHNode = currentCHNode->component;

	for(vector<CH_Node*>::iterator iter = currentCHNode->children.begin();iter!=currentCHNode->children.end();iter++) {
		childrenArea+=(*iter)->traverseAreas();
	}

	tempArea = area;
	if(area<childrenArea)
		tempArea = childrenArea;

	return tempArea;
}

void CH_Node::deleteCh_Node(CH_Node* node) {

	if(node->component!=NULL)
		deleteCh_Node(node->component);

	for(vector<CH_Node*>::iterator iter = node->children.begin();iter!=node->children.end();iter++) {
		deleteCh_Node(*iter);
	}

	delete node;
}

void CH_Node::getClusters(vector<CH_Node*>* clusters) {

	bool thisPushed = false;

	if(component==NULL && children.size()==0) {
		clusters->push_back(this);
		thisPushed = true;
	}

	double area = this->getArea();

	//get areas of the children
	double childrenArea = 0;
	CH_Node* currentCHNode = this;
	while(currentCHNode->component!=NULL)
		currentCHNode = currentCHNode->component;

	for(vector<CH_Node*>::iterator iter = currentCHNode->children.begin();iter!=currentCHNode->children.end();iter++) {
		childrenArea+=(*iter)->tempArea;
	}

	if(area>childrenArea) {
		if(!thisPushed)
			clusters->push_back(this);
	}
	else {
		for(vector<CH_Node*>::iterator iter = currentCHNode->children.begin();iter!=currentCHNode->children.end();iter++) {
			(*iter)->getClusters(clusters);
		}
	}
}


