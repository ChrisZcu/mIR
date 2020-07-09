#ifndef _SKYLINE_H_
#define _SKYLINE_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "filemem.h"
#include "global.h"
#include "virtualRNode.h"
#define MAXPAGEID 49999999

void rtreeRAM(Rtree& rt, unordered_map<long, RtreeNode*>& ramTree);

float minDist(float p1[], float p2[], int dimen);

bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[]);

void GetSkylines(const int dimen, Rtree& a_rtree, std::multimap<long int, VirtualRNode*>& NonResultEntry, std::vector<long int>& PrunedNodes, set<long int>& a_skylines, float* PG[]);

bool isDominateByFocal(const int dimen, const float pt[], Point& focal);

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k);

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k);

void printSkyBand(float* PointSet[], vector<long int>& skyband, int& dim, char* filename);

#endif
