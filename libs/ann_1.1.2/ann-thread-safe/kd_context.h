//----------------------------------------------------------------------
// File:			kd_context.h
// Author:		Andreas Girgensohn
// Description:		Context for search encapsulating previously global variables
//----------------------------------------------------------------------

#ifndef ANN_kd_context_H
#define ANN_kd_context_H

#include <ANN/ANNx.h>					// all ANN includes
#include "pr_queue.h"
#include "pr_queue_k.h"

class ANNkd_context
{
 public:
  int				dim;		// dimension of space (static copy)
  ANNpoint			q;		// query point (static copy)
  double			maxErr;		// max tolerable squared error
  ANNpointArray			pts;		// the points (static copy)
  ANNmin_k			*pointMK;	// set of k closest points
  int				ptsVisited;	// number of points visited
  int				maxPtsVisited;	// maximum number of pts visited

  ANNkd_context(int k) { pointMK = new ANNmin_k(k); }
  ~ANNkd_context() { delete pointMK; }
};

class ANNkd_fr_context : public ANNkd_context
{
 public:
  ANNdist			sqRad;		// squared radius search bound
  int				ptsInRange;	// number of points in the range

  ANNkd_fr_context(int k)  : ANNkd_context(k) {}
};

class ANNkd_pr_context : public ANNkd_context
{
 public:
  double			eps;		// the error bound
  ANNpr_queue			*boxPQ;		// priority queue for boxes

  ANNkd_pr_context(int k, int n_pts) : ANNkd_context(k) {
    boxPQ = new ANNpr_queue(n_pts); // create priority queue for boxes
  }
  ~ANNkd_pr_context() {
    delete boxPQ;		// deallocate priority queue
  }
};
#endif
