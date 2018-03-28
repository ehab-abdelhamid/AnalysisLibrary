/*
 * Clusters.h
 *
 *  Created on: 9 Jan 2018
 *      Author: eabdelha
 */

#ifndef CLUSTERS_H_
#define CLUSTERS_H_

#include<vector>
#include "Cluster.h"

using namespace std;

class Clusters {
private:
	vector<Cluster*> clustersList;

public:

	int getSize() { return clustersList.size(); }

	void addClustrer(Cluster* cluster) {
		cluster->setClusterID(clustersList.size());
		clustersList.push_back(cluster);
	}

	void printElemValues(ostream& out) {
		for(vector<Cluster*>::iterator iter = clustersList.begin();iter!=clustersList.end();iter++) {
			(*iter)->printElemValues(out);
			out<<endl;
		}
	}

	int getClusterID(Object* obj) {
		for(vector<Cluster*>::iterator iter = clustersList.begin();iter!=clustersList.end();iter++) {
			if((*iter)->exists(obj))
				return (*iter)->getClusterID();
		}

		return -1;
	}

	double computePurity(Clusters* otherClusters) {

		double purity = 0;
		int nClusters = otherClusters->getSize();
		int* clusterIDToNumber = new int[nClusters];

		for(vector<Cluster*>::iterator iter = clustersList.begin();iter!=clustersList.end();iter++) {
			Cluster* cluster = (*iter);

			for(int i=0;i<nClusters;i++) {
				clusterIDToNumber[i] = 0;
			}

			int maxValue = 0;
			int maxCluster = -1;
			int count = 0;
			for(vector<Object*>::iterator iter1 = cluster->getElementsBegin();iter1!=cluster->getElementsEnd();iter1++) {
				Object* obj = (*iter1);
				int cid = otherClusters->getClusterID(obj);

				if(cid>=0) {
					clusterIDToNumber[cid]++;
					if(clusterIDToNumber[cid]>maxValue) {
						maxValue = clusterIDToNumber[cid];
						maxCluster = cid;
					}
				}
				count++;
			}


			purity+=(((double)maxValue)/count);
		}

		delete[] clusterIDToNumber;
		purity = purity / clustersList.size();

		return purity;
	}

	~Clusters() {
		for(int i=0;i<clustersList.size();i++) {
			delete clustersList[i];
		}
	}
};

#endif /* CLUSTERS_H_ */
