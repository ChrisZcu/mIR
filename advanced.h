#ifndef ADVANCED_h
#define ADVANCED_H

#include <vector>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/RboxPoints.h>


#include "./rtree/header.h"
#include "./rtree/hypercube.h"
#include "./rtree/rentry.h"
#include "./rtree/rnode.h"
#include "./rtree/rtree.h"
#include "./rtree/filemem.h"
#include "./rtree/tgs.h"
#include "./rtree/global.h"
#include "./rtree/param.h"
#include "./rtree/skyline.h"

#include "./celltree/cellTree.h"
#include "algorithm.h"
#include "util.h"
#include "group.h"

using namespace std;



extern bool union_insert;

extern vector<long int> hsid;

extern vector<vector<float>> user; // user

extern Rtree* rtree; // option (index by R-tree)
struct Cmp{
	bool operator() (const cell *a, const cell *b){
		return a->rank < b->rank;
	}
};

class Advanced:public Algorithm{
public:
      
    Advanced(int dim);
    
    ~Advanced();
   
    double convex_sum;
    double cover_sum;
    double mbr_sum;
    double insert_sum;
    double pre_insert_sum;
    int update_times;
    int update_work;
 
    void solve(int k,int m);

    vector<cell*> look_insert(cell *c, int index);

    void grouping(int k);

    pair<long,int> topkitems(int index,int k);


    vector<int> convex_hull(vector<vector<float>>& points);
   
   int update_group_bo(cell *c,int m);
   vector<cell*> insert(cell *c, int index);

    void cover_check(cell *c, vector<int>& hypers, vector<vector<float>>& vertices, vector<int>& cover, vector<int>& inter, vector<int>& uncover);
    void erase_pruned_groups_users(cell* c, vector<int>& removegroups, unordered_map<int,vector<int>>& covers, unordered_map<int, vector<int>>& uncover);
    void group_cover_check_full(bool& flag, cell *c, int& id, Group& group_copy,vector<int>& hypers, vector<int>& removegroups, vector<vector<float>>& vertices, unordered_map<int,vector<int>>& cover, unordered_map<int,vector<int>>& uncover);
   
    vector<vector<float>> get_vertices(vector<float>& cl, vector<float>& cu);
    
    bool comparepivot(vector<float>& lower, const Point& higher);
    bool comparepivotII(const Point& lower, vector<float>& higher);


    void gen_string(int depth,vector<string>& gen,vector<int>& seq);
    vector<cell*> insert_chull(cell* c, int index);
    vector<cell*> insert_group(cell* c, int index);

   
    void exactly_check(vector<int>& hypers,cell* c,vector<int>& cover, vector<int>& inter, vector<int>& uncover);
     
    vector<float> hyperplane(const Point& option, vector<float>& user);

    priority_queue<cell*,vector<cell*>,Cmp> pq;

private:
    int m_dim;
    int g_cnt;
    int g_prune;
    vector<Group> m_group;  
};






#endif 