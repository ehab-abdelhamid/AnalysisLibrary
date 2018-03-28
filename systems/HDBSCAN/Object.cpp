/*
 * Object.cpp
 *
 *  Created on: 27 Sep 2017
 *      Author: eabdelha
 */

#include <math.h>
#include "Object.h"
#include "Block.h"
#include "BlocksSystem.h"
#include "RandomBlock.h"

long Object::statID = 0;//just modified this on 12-2-2018

void Object::addValue(double a) {
	values.push_back(a);
}

double Object::getValueByIndex(int i) {
	return values[i];
}

void Object::print(ostream& out) {
	cout<<"Values = (";
	for(std::vector<double>::iterator it = values.begin(); it != values.end(); ++it) {
	    out<<*it<<",";
	}
	out<<")"<<endl;
	cout<<"CoreDist = "<<this->coreDist<<endl;
}

void Object::printValues(ostream& out) {
	for(std::vector<double>::iterator it = values.begin(); it != values.end(); ++it) {
	    out<<*it<<",";
	}
}

//double Object::distanceTo(Object* other, double maxDistance) {
//	std::vector<double>::iterator it = values.begin();
//	std::vector<double>::iterator other_it = other->getBeginning();
//
//	double maxDistanceSq = pow(maxDistance, 2);
//
//	double dist = 0;
//
//	for(; it != values.end(); ++it, ++other_it) {
//		dist+=pow(((*it)-(*other_it)), 2);//((*it-*other_it)*(*it-*other_it));//pow((*it-*other_it), 2);
//		if(maxDistance!=-1 && dist>maxDistanceSq) {
//			return maxDistance;
//		}
//	}
//
//	return sqrt(dist);
//}

double Object::distanceTo(Object* other, double maxDistance) {

	double maxDistanceSq = pow(maxDistance, 2);

	double dist = 0;

	for(int i = 0; i<values.size(); i++) {
		dist+=pow((values[i]-other->getValueByIndex(i)), 2);//((*it-*other_it)*(*it-*other_it));//pow((*it-*other_it), 2);
		if(dist>maxDistanceSq && maxDistance!=-1) {
			return maxDistance;
		}
	}

	return sqrt(dist);
}

/**
 * return the max value of core1, core2 and dist(this, other) as the core distance between the two elements
 */
double Object::coreDistanceTo(Object* otherObject, double maxDistance) {
	double dist = this->distanceTo(otherObject, maxDistance);
	if(maxDistance!=-1 && dist>=maxDistance)
		return dist;
	if(dist<this->getCoreDistance())
		dist = this->getCoreDistance();
	if(dist<otherObject->getCoreDistance())
		dist = otherObject->getCoreDistance();
	return dist;
}

double Object::computelength() {

	double length = 0;

	for(std::vector<double>::iterator it = values.begin(); it != values.end(); ++it) {
		length+=((*it)*(*it));
	}

	return sqrt(length);
}

void Object::addValuesFrom(Object* otherO) {

	std::vector<double>::iterator it = this->values.begin();
	for(std::vector<double>::iterator it_other = otherO->getBeginning(); it_other != otherO->getEnd(); ++it_other, ++it) {
		(*it)=(*it)+(*it_other);
	}
}

void Object::setValuesFrom(Object* otherO) {

	for(std::vector<double>::iterator it_other = otherO->getBeginning(); it_other != otherO->getEnd(); ++it_other) {
		this->addValue(*it_other);
	}
}

void Object::divideValuesBy(float v) {
	for(std::vector<double>::iterator it = this->values.begin(); it != values.end(); ++it) {
		(*it)=(*it)/v;
	}
}

void Object::zeroValues(int n) {
	for(int i=0;i<n;i++) {
		values.push_back(0);
	}
}

void Object::calculateCoreDistance(Block* block, BlocksSystem* bs, int nPoints) {
	vector<pair<Object*,double>* >* distances = new vector<pair<Object*,double>* >();
	if(block==NULL)
		block = bs->getBlockOfObject(this);

	block->getCoreDistance(this, nPoints, distances);

	int totalRemaining = nPoints-distances->size();

	double lastDist = distances->back()->second;
	if((lastDist>block->getMarginLength() || distances->size()<nPoints) && bs!=NULL) {
		//iterate over distances, remove ones below margin
		for(vector<pair<Object*,double>* >::iterator iter = distances->begin();iter!=distances->end();) {
			if((*iter)->second<block->getMarginLength()) {
				distances->erase(distances->begin());
				iter = distances->begin();
			}
			else
			{
				iter++;
			}
		}

		int remainingPoints = distances->size()+totalRemaining;
		//go over blocks and try to get ore distances using them
		for(vector<Block*>::iterator iter = bs->getBegin();iter!=bs->getEnd();iter++) {
			if((*iter)==block)
				continue;
			(*iter)->getCoreDistance(this, remainingPoints, distances);
		}
		this->coreDist = distances->back()->second;
	}
	else
	{
		this->coreDist = lastDist;
	}
}

void Object::calculateCoreDistance(RandomBlock* block, int nPoints) {
	vector<pair<Object*,double>* >* distances = new vector<pair<Object*,double>* >();

	block->getCoreDistance(this, nPoints, distances);

	double lastDist = distances->back()->second;
	this->coreDist = lastDist;
}

/**
 * get the K nearest neighbour of this object
 */
pair<Object*, double>* Object::getKNNObject(Block* block, BlocksSystem* bs, int k) {
	vector<pair<Object*,double>* >* distances = new vector<pair<Object*,double>* >();
	if(block==NULL)
		block = bs->getBlockOfObject(this);

	block->getCoreDistance(this, k, distances);

	int totalRemaining = k-distances->size();

	double lastDist = distances->back()->second;
	if((lastDist>block->getMarginLength() || distances->size()<k) && bs!=NULL) {
		//iterate over distances, remove ones below margin
		for(vector<pair<Object*,double>* >::iterator iter = distances->begin();iter!=distances->end();) {
			if((*iter)->second<block->getMarginLength()) {
				distances->erase(distances->begin());
				iter = distances->begin();
			}
			else
			{
				iter++;
			}
		}

		int remainingPoints = distances->size()+totalRemaining;
		//go over blocks and try to get ore distances using them
		for(vector<Block*>::iterator iter = bs->getBegin();iter!=bs->getEnd();iter++) {
			if((*iter)==block)
				continue;
			(*iter)->getCoreDistance(this, remainingPoints, distances);
		}
	}

//	cout<<"Debug: the k nearest object is: ";
//	(distances->back())->first->print(cout);

	return (distances->back());
}

vector<pair<Object*, double>*>* Object::getObjectsWithDistance(Block* block, BlocksSystem* bs, double min, double max) {
	vector<pair<Object*,double>* >* distances = new vector<pair<Object*,double>* >();
	if(block==NULL)
		block = bs->getBlockOfObject(this);

	block->getCoreDistance(this, -1, distances);

	if(min>block->getMarginLength() || min>block->getMarginLength()) {

		//go over blocks and try to get ore distances using them
		for(vector<Block*>::iterator iter = bs->getBegin();iter!=bs->getEnd();iter++) {
			if((*iter)==block)
				continue;

			//create a condition here for pruning blcks with the following condition:
			//(distance(this, (*iter).center)>((*iter).diameter+max) ... think also about what to do for min
			//pruning conditions
			double blocksDistance = (*iter)->getCenterPoint()->distanceTo((*iter)->getCenterPoint());
			double blockRadius = (*iter)->getRadius();
			if(max<blocksDistance-blockRadius || min>blocksDistance+blockRadius)
				cout<<"*";
			else
				(*iter)->getCoreDistance(this, -1, distances);
		}
	}

	//remove objects with higher or lower distances than needed
	vector<pair<Object*,double>* >::iterator iter = distances->begin();
	for(;iter!=distances->end();) {
		pair<Object*,double>* p_o = (*iter);
//		p_o->first->print(cout);
//		cout<<" --- ";
//		this->print(cout);
//		cout<<" .distance = ";
//		cout<<p_o->second<<","<<min<<","<<max<<endl;
		if(p_o->second<min || p_o->second>max)
			iter = distances->erase(iter);
		else
			iter++;
	}

	return distances;
}

/**
 * get the kth shortest reachability distance
 * if core distance is not calculated, then calculate it using coreK
 */
pair<Object*, double>* Object::getKShortestRechability(Block* block, BlocksSystem* bs, int coreK, int k) {
	if(block==NULL)
		block = bs->getBlockOfObject(this);

	if(coreDist==-1)
		calculateCoreDistance(block, bs, coreK);

	pair<Object*, double>* po = getKNNObject(block, bs, k);

	double otherCoreDistance = po->first->getCoreDistance();
	if(otherCoreDistance==-1) {
		po->first->calculateCoreDistance(NULL, bs, coreK);
		otherCoreDistance = po->first->getCoreDistance();
	}

	//get the max of this->coreDist, other->coreDist and distance
	double maxDist = coreDist;

	if(otherCoreDistance>maxDist)
		maxDist = otherCoreDistance;

	if(po->second>maxDist)
		maxDist = po->second;

	return new pair<Object*, double>(po->first, maxDist);
}
