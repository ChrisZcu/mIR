#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>

#include "./celltree/cellTree.h"


class Algorithm{

public:
    Algorithm();
    vector<cell*> m_final;
    cellTree* m_celltree;
    virtual void solve(int k,int m) = 0;
protected:
    
};

#endif