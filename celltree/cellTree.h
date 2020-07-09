/* ----------------------------------------------------------------------------
This header file includes class Advanced RegionTree declaration.
It is data structure to maintain the cells in query plane with hyperplanes
---------------------------------------------------------------------------- */

#ifndef _CELL_TREE_H_
#define _CELL_TREE_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>


#include "lp_lib.h"
#include "../group.h"
#include "../rtree/rtree.h"
#include "../rtree/rentry.h"
#include "../rtree/rnode.h"
#include "../rtree/filemem.h"


using namespace std;
extern vector<vector<float>> HalfSpaces;
extern map<long int, long int> RecordIDtoHalfPlaneID;
extern long int totalSpaceCost;
extern unordered_map<long int, RtreeNode*> ramTree;
extern int objCnt;
extern int node_size;
extern int m;
extern int num_cell;
extern int num_lp;

struct rBound
{
	int minimum;
	int maximum;

	rBound()
	{
		minimum = 0;
		maximum = objCnt - 1;
	}
};

struct cell
{
	unordered_map<long int, bool> IntersectHP;
	unordered_set<long int> AboveHP;
	unordered_set<long int> BelowHP;

	vector<int> G;
	rBound lu;
	int rank;
	bool isPruned;
	int lower;
	int upper;
	int depth;
	vector<float> cu,cl;
	cell* left;
	cell* right;
	cell()
	{	
		depth = 0;
		num_cell++;
		left = NULL;
		right = NULL;
		isPruned = false;
		rank = 0;
		lower = 0;
		upper = 0;
	}

	cell(cell* copy)
	{
		
			depth = copy->depth + 1;
			num_cell++;
			left = NULL;
			right = NULL;
			IntersectHP = copy->IntersectHP;
			AboveHP = copy->AboveHP;
			BelowHP = copy->BelowHP;
			rank = copy->rank;
			isPruned = copy->isPruned;
			G = copy->G; // important
			lower = copy->lower;
			upper = copy->upper;
	}

	void appendto(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
			node->IntersectHP[iter->first] = iter->second;
		}

		for (auto iter = AboveHP.begin(); iter != AboveHP.end(); iter++)
		{
			node->AboveHP.insert(*iter);
		}

		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void copyleaf(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
				node->IntersectHP[iter->first] = iter->second;
			
		}
		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		for (auto iter = AboveHP.begin(); iter != AboveHP.end(); iter++)
		{
			node->AboveHP.insert(*iter);
		}
		node->G = G;
		node->lower = lower;
		node->upper = upper;
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void release()
	{
		IntersectHP.clear();
		unordered_map<long int, bool>().swap(IntersectHP);
		AboveHP.clear();
		unordered_set<long int>().swap(AboveHP);
		BelowHP.clear();
		unordered_set<long int>().swap(BelowHP);
	}

};

struct gNode
{
	long int rID;
	unordered_set<long int> rDominator;
	gNode(long int recordID) :rID(recordID){};
};

struct cellCompare
{
	bool operator()(const cell* a, const cell* b) const
	{
		return a->rank < b->rank;
	}
};


class cellTree
{
public:
	cellTree();
	~cellTree();

	int chull_size;  
	int hyper_size; // update the rank if a cell covered by intersect of convex hull
	
	map<cell*,int> inserted;

	void releaseCell(cell* node);

	void insert(cell *c, vector<long int> &hps, int mink, vector<cell*>& leaves);//, vector<cell*>& finalResult);
	bool isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator);
	bool isCellFeasible(unordered_map<long int, bool>& touchhs);
	void lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs);
	void addHP(lprec* model, long int hpid, bool sideindicator);

	void updateRank(cell* node, const int mink);
	void updateUpper(cell* node, const int mink);
	void inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs);

	void collectLeaf(cell *c, vector<cell*>& leaves, const int& mink);
	void dfsTraversal(cell* node, vector<cell*>& leaves);


	void opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult);
	void opt_inserthp(long int & hpid, const int mink, cell* node, cell* all);
	int isDominatorInserted(unordered_set<long int>& rdominators, cell* all);


	void updateCellTree(cell* node);


	int findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu);

	cell* root;
	unordered_map<long int, gNode*> dagNode;
};



#endif