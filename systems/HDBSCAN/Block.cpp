/*
 * Block.cpp
 *
 *  Created on: 3 Oct 2017
 *      Author: eabdelha
 */

#include<cstdlib>
#include "Block.h"
#include "Utils.h"
#include "BlocksSystem.h""


Block::Block() {
	marginLength = 100;
	left = NULL;
	right = NULL;
}

void Block::addElement(Object* element, bool coreElement) {
	if(coreElement) {
		elements.insert(std::make_pair<int, Object*>(element->getID(), element));
		element->setBlock(this);
	}
	else
		borderElems.insert(std::make_pair<int, Object*>(element->getID(), element));
}

void Block::addElemets(vector<Object*> elements, bool coreElement) {
	for(vector<Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
		addElement(*iter, coreElement);
	}
}

Object* Block::getElementByID(int id) {
	map<int, Object*>::iterator iter = elements.find(id);
	if(iter==elements.end())
		return NULL;
	else
		return iter->second;
}

/**
 * pick a random element from the list of elements
 */
Object* Block::pickRandomElement() {
	int r_index = rand() % elements.size();
	map<int, Object*>::iterator iter = elements.begin();
	for(int i=0;i<r_index;i++)
		iter++;
	return iter->second;
	//return elements.at(r_index);
}

/**
 * select the node which is the farthest to the given point
 */
Object* Block::getFarthestPoint(Object* obj) {

	Object* r_obj = obj;
	double max_distance = 0;

	for(std::map<int, Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		Object* curr_obj = it->second;
		double curr_distance = obj->distanceTo(curr_obj);
		if(curr_distance>max_distance) {
			max_distance = curr_distance;
			r_obj = curr_obj;
		}
	}

	return r_obj;
}

void Block::createProjectionLine() {

	Object* r = pickRandomElement();	//pick a random point
	Object* p1 = getFarthestPoint(r);	//select the farthest point from the random one
	Object* p2 = getFarthestPoint(p1);	//select the farthest point from p1

	if(p1==p2) {
		cout<<"p1==p2!!!"<<endl;
		cout<<"size = "<<this->elements.size()<<endl;
		p2 = getFarthestPoint(p1);
	}

	projectionLine = new pair<Object*, Object*>(p1, p2);
}

/**
 * project the given point on the projection line, return the resulting point
 */
Object* Block::projectOnProjectionline(Object* obj) {
	if(this->projectionLine==NULL)
		createProjectionLine();

	return getProjection(obj, projectionLine);
}

/**
 * do the projection given the point and projection line
 */
Object* Block::getProjection(Object* obj, std::pair<Object*, Object*>* projectionLine) {

	// get dot product of e1, e2
	Object* v1 = projectionLine->first;
	Object* v2 = projectionLine->second;

	if(obj==v1 || obj==v2)
	{
		Object* p = new Object();
		std::vector<double>::iterator it = obj->getBeginning();

		for(;it!=obj->getEnd();++it) {
			p->addValue(*it);
		}
		return p;
	}

	Object* e1 = subtract(v2, v1);
	Object* e2 = subtract(obj, v1);

	double valDp = dot(e1, e2);

	// get length of vectors
	double lenLineE1 = e1->computelength();
	double lenLineE2 = e2->computelength();

	double cos = valDp / (lenLineE1 * lenLineE2);
	// length of v1P'
	double projLenOfLine = cos * lenLineE2;

	Object* p = new Object();
	std::vector<double>::iterator it_v1 = v1->getBeginning();
	std::vector<double>::iterator it_e1 = e1->getBeginning();

	for(;it_v1!=v1->getEnd();++it_v1,++it_e1) {
		p->addValue((*it_v1)+(projLenOfLine * (*it_e1))/lenLineE1);
	}

	return p;
}

/**
 * divide the dataset into two overlapping parts; left and right
 */
void Block::partition(vector<Block*>* smallestBlocks) {

	//stop condition: you can modify it
	if(this->elements.size()<=100)
	{
		cout<<"No more partitions, size is "<<elements.size()<<" which is small."<<endl;
		if(smallestBlocks!=NULL)
			smallestBlocks->push_back(this);
		return;
	}

	right = new Block();
	left = new Block();

	createProjectionLine();

	//project objects on the projection line
	vector<Object*> projectedElements;
	for(std::map<int, Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		Object* curr_obj = it->second;
		Object* preojectedObject = projectOnProjectionline(curr_obj);
		projectedElements.push_back(preojectedObject);
	}

	//get the object on the middle of the line
	Object* middle = new Object();
	Object* v1 = projectionLine->first;
	Object* v2 = projectionLine->second;
	std::vector<double>::iterator it2 = v2->getBeginning();
	for(std::vector<double>::iterator it1 = v1->getBeginning(); it1 != v1->getEnd(); ++it1, ++it2) {
		middle->addValue((*it1+*it2)/2.0);
	}

	//compute the extended boundary points
	double d_ = v1->distanceTo(v2);
	double d = marginLength;
	//get the point on the right (d=+ve) and the left point (-ve)
	Object* rightB = new Object();
	Object* leftB = new Object();
	std::vector<double>::iterator it1 = v1->getBeginning();
	it2 = v2->getBeginning();
	for(std::vector<double>::iterator it = middle->getBeginning(); it != middle->getEnd(); ++it,++it1,++it2) {
		double newValue = (*it)+d;//+(d/d_)*((*it2)-(*it1));
		rightB->addValue(newValue);
		newValue = (*it)-d;//+(-1*d/d_)*((*it2)-(*it1));
		leftB->addValue(newValue);
	}

	//partition objects based on the middle node
	//first, select a dimension on which we compare the objects.
	//The only requirement for a dimension is to have a difference between v1 and v2 in this dimension
	int dIndex = 0;//the selected dimension index
	while(v1->getValueByIndex(dIndex)==v2->getValueByIndex(dIndex))
		dIndex++;
//	bool positive = true;
//	if((v2->getValueByIndex(dIndex)-v1->getValueByIndex(dIndex))<0)
//		positive = false;

	int counter = 0;
	std::map<int, Object*>::iterator it = elements.begin();
	for(std::vector<Object*>::iterator projectedIt = projectedElements.begin(); projectedIt != projectedElements.end(); ++projectedIt,++it) {
		counter++;
		Object* curr_obj = *projectedIt;

		//adding core elements
		if(curr_obj->getValueByIndex(dIndex)>middle->getValueByIndex(dIndex))
			right->addElement(it->second, true);
		else
			left->addElement(it->second, true);

		//adding boundary elements
		if(curr_obj->getValueByIndex(dIndex)>middle->getValueByIndex(dIndex) && curr_obj->getValueByIndex(dIndex)<=rightB->getValueByIndex(dIndex))
			left->addElement(it->second, false);
		if(curr_obj->getValueByIndex(dIndex)<=middle->getValueByIndex(dIndex) && curr_obj->getValueByIndex(dIndex)>=leftB->getValueByIndex(dIndex))
			right->addElement(it->second, false);

		/*if(positive) {
			if(curr_obj->getValueByIndex(dIndex)<=rightB->getValueByIndex(dIndex))
				right->addElement(it->second);
			if(curr_obj->getValueByIndex(dIndex)>=leftB->getValueByIndex(dIndex))
				left->addElement(it->second);
		} else {
			if(curr_obj->getValueByIndex(dIndex)>=rightB->getValueByIndex(dIndex))
				right->addElement(it->second);
			if(curr_obj->getValueByIndex(dIndex)<=leftB->getValueByIndex(dIndex))
				left->addElement(it->second);
		}*/
	}

	cout<<"Block of size "<<elements.size()<<" is partitioned into two partitions with sizes: "<<right->getNumberOfElements()<<", "<<left->getNumberOfElements()<<endl;

	//then recursively partition right and left
	right->partition(smallestBlocks);
	left->partition(smallestBlocks);
}

void Block::calcCentreAndRadius() {

	//comput the center point
	centerPoint = new Object();

	if(elements.size()==0)
	{
		return;
	}

	std::map<int, Object*>::iterator it = elements.begin();
	centerPoint->setValuesFrom(it->second);
	++it;

	for(; it != elements.end(); ++it) {
		centerPoint->addValuesFrom(it->second);
	}

	centerPoint->divideValuesBy(elements.size());

	//compute the radius (the maximum distance from center point to points)
	this->radius = 0;
	for(std::map<int, Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		double distance = centerPoint->distanceTo(it->second);
		if(distance>radius)
			radius = distance;
	}
}

void Block::computeCoreDistances(BlocksSystem* bs) {
	for(std::map<int, Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		Object* object = it->second;
		object->calculateCoreDistance(this, bs, difaultNPoints);
	}
}

void Block::getCoreDistance(Object* obj, int knn, vector<pair<Object*,double>* >* sortedList) {

	map<int, Object*>* list = &elements;

	//go over the two lists; elements and borderElems
	while(true) {

		for(std::map<int, Object*>::iterator it = list->begin(); it != list->end(); ++it) {
			Object* otherObject = it->second;
			if(obj==otherObject)
				continue;

			double distance = obj->distanceTo(otherObject);
			//put it in order
			int count = 0;
			bool insert = true;
			std::vector<pair<Object*,double>*>::iterator its = sortedList->begin();
			for(; its != sortedList->end(); ++its) {
				pair<Object*, double>* o = *its;
				if(o->first==otherObject) {
					insert = false;
					break;
				}
				if(o->second>distance)
					break;
				count++;
				if(knn!=-1 && count>knn)
					break;
			}

			if(insert)
				sortedList->insert(its, new pair<Object*, double>(otherObject, distance));

			if(knn!=-1 && sortedList->size()>knn)
			{
				delete sortedList->back();
				sortedList->pop_back();
			}
		}

		if(list==(&elements))
			list = &borderElems;
		else
			break;

	}
}

set<Object* >* Block::getObjectsSet() {
	set<Object*>* s = new set<Object*>();

	for(std::map<int, Object*>::iterator it = elements.begin(); it != elements.end(); ++it) {
		Object* o = it->second;
		s->insert(o);
	}

	return s;
}

/**
 * THINK!
 * classifies an object whether it should go to the right or left partition
 * The object is added either to the left or right block
 * returns 1 of right, 2 of left
 */
int Block::classify(Object*) {

}

void Block::print(ostream& out, bool printElems) {
	cout<<"#Core Objects = "<<elements.size()<<", #Border Objects = "<<borderElems.size()<<endl;
	if(centerPoint!=NULL) {cout<<"Center = ";centerPoint->print(cout);}
	cout<<"Radius = "<<this->radius<<endl;

	if(printElems) {
		out<<"Elements:"<<endl;
		for(map<int, Object*>::iterator iter = elements.begin();iter!=elements.end();iter++) {
			iter->second->print(out);
		}
	}
}

Block::~Block() {
	if(projectionLine!=NULL) delete projectionLine;
}
