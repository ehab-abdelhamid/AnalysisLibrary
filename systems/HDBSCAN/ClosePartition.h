#ifndef ClosePartition_H_
#define ClosePartition_H_

#include<vector>
#include"Object.h"

using namespace std;

class PossibleEdge;

class ClosePartition {
public:
	static vector<ClosePartition*> generatePartitions(vector<Object*>, double r);
	bool isWithinDist(Object* obj);
	double getDistanceToCenter(Object* obj);
	bool addObject(Object* obj);
	void print(ostream& out);
	void setR(double r) { this->r = r; }
	long getR() { return r; }
	vector<PossibleEdge*> getElementsWithDistBound(Object* obj, double min, double max);

private:
	long r = -1;
	Object* centroid;
	vector<Object*> objects;
};

#endif /* ClosePartition_H_ */
