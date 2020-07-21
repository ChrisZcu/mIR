#ifndef BASELINE_H
#define BASELINE_H

#include <vector>
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
#include "algorithm.h"
#include "product.h"
#include "./celltree/cellTree.h"

extern vector<long int> hsid;
extern Rtree* rtree;
extern vector<vector<float>> HalfSpaces;


class BaseLine:public Algorithm{
public:
    vector<cell*> m_final;
    BaseLine(//const unordered_map<long, RtreeNode *>& ramTree
            const std::vector<User>& m_user,int dim);

    ~BaseLine();
    Product solve(int k,double ratio);
    //first we get all top-k hyperplanes of all users
    void gettopk(int k);
    //insert planes into product space
    void insert(int m);
    // minize cost in all avaliable space
    void optimize();
private:
    std::vector<User> m_user;
    Rtree* m_rtree; 
    int m_dim;
    vector<pair<double,Product>> m_result;
    cellTree* m_celltree;
  //  unordered_map<long, RtreeNode *> m_ramTree;

};


#endif