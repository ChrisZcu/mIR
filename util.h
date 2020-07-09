#ifndef UTIL_H
#define UTIL_H
#include <cstring>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <queue>
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
#include "celltree/cellTree.h"


extern vector<vector<float>> HalfSpaces;
//read parameter from command line
char* param(int argc,char **argv,char *target,char *def);

void load_user(vector<vector<float>>& user,int dim,const char* path,int n);
    
void testtopk(const char* datafile,int k,int dim,int psize,const vector<vector<float>>& user);

void display_cell(std::vector<cell*>& cells);

double user_score(const vector<float>& user,const Point& p);

double user_score(const vector<float>& user,const vector<float>& product);

void visit_hs(int id);
#endif