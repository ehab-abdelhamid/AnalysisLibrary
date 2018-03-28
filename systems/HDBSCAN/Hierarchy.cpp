/*
 * Hierarchy.cpp
 *
 *  Created on: 9 Jan 2018
 *      Author: eabdelha
 */

#include "Hierarchy.h"

struct Qu_mcp {
public:
	H_Node* ths;
	int minSize;
	CH_Node* ch;
	H_Node* parent;
	Qu_mcp(H_Node* t, int ms, CH_Node* ch, H_Node* p) {ths = t;minSize = ms;this->ch = ch;parent = p;}
};

int H_Node::getNumItems() {
	return numItems;
}

set<NodeX*>* H_Node::getItems() {

	set<NodeX*>* items = new set<NodeX*>();

	//collect all nodes in the tree in
	vector<std::pair<H_Node*, H_Node*>*> toBeProcessed;
	toBeProcessed.push_back(new std::pair<H_Node*, H_Node*>(this, NULL));

	while(toBeProcessed.size()>0) {
		std::pair<H_Node*, H_Node*>* currentPair = *toBeProcessed.begin();
		H_Node* thisNode = currentPair->first;
		//H_Node* parent = currentPair->second;
		toBeProcessed.erase(toBeProcessed.begin());
		delete currentPair;

		if(thisNode->nodeValue!=NULL)
			items->insert(thisNode->nodeValue);

		if(thisNode->elements.size()==0 || (thisNode->elements.size()==1 && (*thisNode->elements.begin())->hnode==thisNode->parent)) {
		}
		else {
			for(vector<HNode_Distance*>::iterator iter = thisNode->elements.begin();iter!=thisNode->elements.end();iter++) {
				//do not iterate over the parent
				if((*iter)->hnode==thisNode->parent)
					continue;

				toBeProcessed.push_back(new std::pair<H_Node*, H_Node*>((*iter)->hnode, thisNode));
			}
		}
	}

	return items;
}

void H_Node::setNumItems(H_Node* p) {

	//collect all nodes in the tree in
	vector<std::pair<H_Node*, H_Node*>*> toBeProcessed;
	vector<std::pair<H_Node*, H_Node*>*> total;
	toBeProcessed.push_back(new std::pair<H_Node*, H_Node*>(this, p));
	this->parent = p;

	while(toBeProcessed.size()>0) {
//		cout<<"toBeProcessed.size() = "<<toBeProcessed.size()<<endl<<flush;
		std::pair<H_Node*, H_Node*>* currentPair = *toBeProcessed.begin();
		H_Node* thisNode = currentPair->first;
		H_Node* parent = currentPair->second;
		toBeProcessed.erase(toBeProcessed.begin());
		delete currentPair;

		if(thisNode->elements.size()==0 || (thisNode->elements.size()==1 && (*thisNode->elements.begin())->hnode==parent)) {
		}
		else {
			for(vector<HNode_Distance*>::iterator iter = thisNode->elements.begin();iter!=thisNode->elements.end();iter++) {
				//do not iterate over the parent
				if((*iter)->hnode==parent)
					continue;

				toBeProcessed.push_back(new std::pair<H_Node*, H_Node*>((*iter)->hnode, thisNode));
				(*iter)->hnode->parent = thisNode;
				total.push_back(new std::pair<H_Node*, H_Node*>((*iter)->hnode, thisNode));
			}
		}
	}

	while(total.size()>0) {
		std::pair<H_Node*, H_Node*>* currentPair = total.back();
		H_Node* thisNode = currentPair->first;
		H_Node* parent = currentPair->second;
		total.pop_back();
		delete currentPair;

		if(thisNode->elements.size()==0 || (thisNode->elements.size()==1 && (*thisNode->elements.begin())->hnode==parent)) {
			if(thisNode->nodeValue==NULL) {
				thisNode->numItems = 0;
			}
			else {
				thisNode->numItems = 1;
			}
		}
		else {
			thisNode->numItems = 0;
			for(vector<HNode_Distance*>::iterator iter = thisNode->elements.begin();iter!=thisNode->elements.end();iter++) {
				//do not iterate over the parent
				if((*iter)->hnode==parent)
					continue;

				thisNode->numItems+=(*iter)->hnode->numItems;
			}
		}

//		cout<<thisNode<<", "<<thisNode->numItems<<endl<<flush;
	}
}

/*void H_Node::setNumItems(H_Node* parent) {
	if(elements.size()==0 || (elements.size()==1 && (*elements.begin())->hnode==parent)) {
		if(nodeValue==NULL) {
			numItems = 0;
		}
		else {
			numItems = 1;
		}
	}
	else {
		numItems = 0;
		for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
			//do not iterate over the parent
			if((*iter)->hnode==parent)
				continue;
			(*iter)->hnode->setNumItems(this);
			numItems+=(*iter)->hnode->getNumItems();
			items.insert((*iter)->hnode->getItems()->begin(), (*iter)->hnode->getItems()->end());
		}
	}
}*/

void H_Node::addNode(H_Node* hn, double distance) {
	elements.push_back(new HNode_Distance(hn, distance));
}

void H_Node::setNode(NodeX* n) {
	nodeValue = n;
//	items.insert(n);
}

int H_Node::getNodeID() {
	if(nodeValue==NULL)
		return -1;
	else {
		return nodeValue->getID();
	}
}

void H_Node::condense(int ms, CH_Node* c, H_Node* p) {

	vector<Qu_mcp*> toBeProcessed;
	Qu_mcp* mcp = new Qu_mcp(this, ms, c, p);
	toBeProcessed.push_back(mcp);

	while(toBeProcessed.size()>0) {

//		if(toBeProcessed.size()%1000000==0)
//			cout<<"toBeProcessed.size() = "<<toBeProcessed.size()<<endl<<flush;
		mcp = *(toBeProcessed.begin());
		H_Node* thisNode = mcp->ths;
		int minSize = mcp->minSize;
		CH_Node* ch = mcp->ch;
		H_Node* parent = mcp->parent;
		delete mcp;//this is causing the problem!
		toBeProcessed.erase(toBeProcessed.begin());

		ch->setDataNode(thisNode);
		if(thisNode->elements.size()==0 || (thisNode->elements.size()==1 && (*thisNode->elements.begin())->hnode==parent)) {
			continue;
		}

		double endValue = -1;
		vector<HNode_Distance*> forRecurion;

		//get the next edge to remove
		for(vector<HNode_Distance*>::iterator iter = thisNode->elements.begin();iter!=thisNode->elements.end();iter++) {

			double tempDist = 1.0/(*iter)->distance;

			if((*iter)->hnode==parent)
				continue;

			if(endValue!=-1 && endValue!=tempDist) {
				cout<<"Errorrrrrrrrrrrrrrrrrrrrrrrr (126172)"<<endl;
			}
			endValue = tempDist;

			if((*iter)->hnode->getNumItems()<minSize) {
				continue;
			}

			forRecurion.push_back((*iter));
		}

		ch->setEnd(endValue);

		if(forRecurion.size()==1) {
			//discard the small one, recursively condense the remaining
			CH_Node* nch = new CH_Node();
			nch->setStart(endValue);
			nch->setWidth((*forRecurion.begin())->hnode->getNumItems());
			ch->setComponent(nch);

			Qu_mcp* nmcp = new Qu_mcp((*forRecurion.begin())->hnode, minSize, nch, thisNode);
			toBeProcessed.push_back(nmcp);

			//(*forRecurion.begin())->hnode->condense(minSize, nch, this);
		} else if(forRecurion.size()==2) {
			HNode_Distance* elem1 = (*forRecurion.begin());
			HNode_Distance* elem2 = (*(forRecurion.begin()+1));
			int num1 = elem1->hnode->getNumItems();
			int num2 = elem2->hnode->getNumItems();
			//create two children, one for the new component and one for the remaining
			//the first child; the one following th enew link
			CH_Node* nch1 = new CH_Node();
			nch1->setStart(endValue);
			nch1->setWidth(num1);
			ch->addChild(nch1);

			Qu_mcp* nmcp = new Qu_mcp(elem1->hnode, minSize, nch1, thisNode);
			toBeProcessed.push_back(nmcp);
			//elem1->hnode->condense(minSize, nch1, this);

			//the second child; the remaining
			CH_Node* nch2 = new CH_Node();
			nch2->setStart(endValue);
			nch2->setWidth(num2);
			ch->addChild(nch2);

			nmcp = new Qu_mcp(elem2->hnode, minSize, nch2, thisNode);
			toBeProcessed.push_back(nmcp);
			//elem2->hnode->condense(minSize, nch2, this);
		} else if(forRecurion.size()>2) {
			cout<<"Errorrrrrrrrrrrrrrrrrrrrrrrr (6843854)"<<endl;
		}
	}
}

/**
 * Old recursive condense function. It is correct
 */
/*void H_Node::condense(int minSize, CH_Node* ch, H_Node* parent) {

	ch->setDataNode(this);

	if(elements.size()==0 || (elements.size()==1 && (*elements.begin())->hnode==parent)) {
		return;
	}

	double endValue = -1;
	vector<HNode_Distance*> forRecurion;

	//get the next edge to remove
	for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {

		double tempDist = 1.0/(*iter)->distance;

		if((*iter)->hnode==parent)
			continue;

		if(endValue!=-1 && endValue!=tempDist) {
			cout<<"Errorrrrrrrrrrrrrrrrrrrrrrrr (126172)"<<endl;
		}
		endValue = tempDist;

		if((*iter)->hnode->getNumItems()<minSize) {
			continue;
		}

		forRecurion.push_back((*iter));
	}

	ch->setEnd(endValue);

	if(forRecurion.size()==1) {
		//discard the small one, recursively condense the remaining
		CH_Node* nch = new CH_Node();
		nch->setStart(endValue);
		nch->setWidth((*forRecurion.begin())->hnode->getNumItems());
		ch->setComponent(nch);
		(*forRecurion.begin())->hnode->condense(minSize, nch, this);
	} else if(forRecurion.size()==2) {
		HNode_Distance* elem1 = (*forRecurion.begin());
		HNode_Distance* elem2 = (*(forRecurion.begin()+1));
		int num1 = elem1->hnode->getNumItems();
		int num2 = elem2->hnode->getNumItems();
		//create two children, one for the new component and one for the remaining
		//the first child; the one following th enew link
		CH_Node* nch1 = new CH_Node();
		nch1->setStart(endValue);
		nch1->setWidth(num1);
		ch->addChild(nch1);
		elem1->hnode->condense(minSize, nch1, this);

		//the second child; the remaining
		CH_Node* nch2 = new CH_Node();
		nch2->setStart(endValue);
		nch2->setWidth(num2);
		ch->addChild(nch2);
		elem2->hnode->condense(minSize, nch2, this);
	} else if(forRecurion.size()>2) {
		cout<<"Errorrrrrrrrrrrrrrrrrrrrrrrr (6843854)"<<endl;
	}
}*/

/*void H_Node::condense(int minSize, CH_Node* ch, H_Node* parent = NULL) {

	if(elements.size()==0 || (elements.size()==1 && (*elements.begin())->hnode==parent))
		return;

	//flag to indicate that there is an element with minimum size
	bool minSizeExists = false;

	//get the next edge to remove
	vector<HNode_Distance*>::iterator shortestEdge = elements.end();
	double shortestDist = 1000000;
	for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {

		double tempDist = 1.0/(*iter)->distance;

		if((*iter)->hnode->getNumItems()<minSize) {
			minSizeExists = true;
			if(shortestEdge==elements.end())
				ch->setEnd(tempDist);
			continue;
		}

		if((*iter)->hnode==parent)
			continue;

		if(tempDist<shortestDist) {
			shortestDist = tempDist;
			shortestEdge = iter;
		}
	}

	if(shortestEdge==elements.end())
		return;

	HNode_Distance* toRemoveCHNode = (*shortestEdge);
	elements.erase(shortestEdge);
	ch->setEnd(shortestDist);

	if(minSizeExists) {
		//discard the small one, recursively condense the remaining
		CH_Node* nch = new CH_Node();
		nch->setStart(shortestDist);
		nch->setWidth(toRemoveCHNode->hnode->getNumItems());
		ch->setComponent(nch);
		toRemoveCHNode->hnode->condense(minSize, nch, this);
	} else {
		int num1 = toRemoveCHNode->hnode->getNumItems();
		int num2 = ch->getWidth()-toRemoveCHNode->hnode->getNumItems();
		//create two children, one for the new component and one for the remaining
		//the first child; the one following th enew link
		CH_Node* nch1 = new CH_Node();
		nch1->setStart(shortestDist);
		nch1->setWidth(num1);
		ch->addChild(nch1);
		toRemoveCHNode->hnode->condense(minSize, nch1, this);

		//the second child; the remaining
		CH_Node* nch2 = new CH_Node();
		nch2->setStart(shortestDist);
		nch2->setWidth(num2);
		ch->addChild(nch2);
		condense(minSize, nch2, this);
	}
}*/

/*void condense(int minSize, CH_Node* ch, H_Node* parent = NULL) {
	double tempDist;
	for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();) {
		int size = (*iter)->hnode->getNumItems();
		tempDist = 1.0/(*iter)->distance;
		if(size<minSize) {
			iter = elements.erase(iter);
		} else {
			iter++;
		}
	}

	if(elements.size()==0) {
		ch->setEnd(tempDist);
	} else if(elements.size()==1) {
		if((*elements.begin())->hnode==parent) {
			ch->setEnd(tempDist);
		}
		else {
			CH_Node* nch = new CH_Node();
			double currentDistance = 1.0/(*elements.begin())->distance;
			ch->setEnd(currentDistance);
			nch->setStart(currentDistance);
			nch->setWidth((*elements.begin())->hnode->getNumItems());
			ch->setComponent(nch);
			(*elements.begin())->hnode->condense(minSize, nch, this);
		}
	} else {
		for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
			//do not go back to the parent
			if((*iter)->hnode==parent)
				continue;

			CH_Node* nch = new CH_Node();
			double currentDistance = 1.0/(*iter)->distance;
			ch->setEnd(currentDistance);
			nch->setStart(currentDistance);
			nch->setWidth((*iter)->hnode->getNumItems());
			ch->addChild(nch);
			(*iter)->hnode->condense(minSize, nch, this);
		}
	}
}*/

void H_Node::print(H_Node* parent) {

	cout<<"This function is commented, get back to it if needed!"<<endl;
	/*

	if(elements.size()==0 || (elements.size()==1 && (*elements.begin())->hnode==parent)) {
		if(nodeValue==NULL) {
			cout<<"-1"<<endl;
		}
		else {
			cout<<nodeValue->getID()<<endl;
		}
	}
	else {
		cout<<"{";
		for(set<NodeX*>::iterator iter = items.begin();iter!=items.end();iter++) {
			cout<<(*iter)->getID()<<",";
		}
		cout<<"}";
		for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
			//do not go back to the parent
			if((*iter)->hnode==parent)
				continue;
			cout<<"| ["<<(*iter)->distance<<"]";
			(*iter)->hnode->print(this);
		}
	}*/
}

void H_Node::printAsList(H_Node* parent) {
	set<NodeX*>* items = getItems();

	for(set<NodeX*>::iterator iter=items->begin();iter!=items->end();iter++) {
		NodeX* node = (*iter);
		cout<<node->getID()<<",";
		node->getObject()->printValues(cout);
		cout<<"*";
	}

	delete items;
}

void H_Node::deleteDescendents() {

	//collect all nodes in the tree in
	vector<H_Node*> toBeProcessed;
	toBeProcessed.push_back(this);

	while(toBeProcessed.size()>0) {
		H_Node* thisNode = *toBeProcessed.begin();
		toBeProcessed.erase(toBeProcessed.begin());

		if(thisNode->elements.size()==0 || (thisNode->elements.size()==1 && (*thisNode->elements.begin())->hnode==thisNode->parent)) {
		}
		else {
			for(vector<HNode_Distance*>::iterator iter = thisNode->elements.begin();iter!=thisNode->elements.end();iter++) {
				//do not iterate over the parent
				if((*iter)->hnode==thisNode->parent)
					continue;

				toBeProcessed.push_back((*iter)->hnode);
			}
		}

		if(thisNode!=this)
			delete thisNode;
	}
}

H_Node::~H_Node() {

	for(vector<HNode_Distance*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		delete (*iter);
	}
}
