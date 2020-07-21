#include <iostream>
#include <cstring>
#include <fstream>
#include <queue>
#include "util.h"
#include "advanced.h"
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
#include "baseline.h"

vector<vector<float>> user;
long int totalSpaceCost;
int node_size = 0;
vector<vector<float>> HalfSpaces;
map<long int, long int> RecordIDtoHalfPlaneID;

vector<long int> hsid; 
int objCnt = 0; /// object count, init with 0
float** PointSet ;//= new float*[MAXPTS + 1];
RtreeNodeEntry** p;// = new RtreeNodeEntry*[MAXPTS];
unordered_map<long, RtreeNode *> ramTree; // load Rtree to main-memory
double _totalSpaceCost = 0.0; // space cost (MB)
FileMemory* mem;
Rtree* rtree;
int num_cell = 0;
int num_lp = 0;
int m; //  target user number

int objcnt = 0;



void build_rtree(const char* datafile, const char* indexfile,int dim,FileMemory* mem){
    
    cout << "Load data points from file" << endl;
	PointSet = new float*[MAXPTS + 1];
	p = new RtreeNodeEntry*[MAXPTS];
	fstream fpdata;
	fpdata.open(datafile, ios::in);
	while (true)
	{
		int id;
		float* cl = new float[dim];
		float* cu = new float[dim];
		fpdata >> id;
		if (fpdata.eof())
			break;

		PointSet[objCnt + 1] = new float[2 * dim];

		for (int d = 0; d < dim; d++)
		{
			fpdata >> cl[d];
        	PointSet[objCnt + 1][d] = cl[d];
		}

		for (int d = 0; d < dim; d++)
		{
			fpdata >> cu[d];
        	PointSet[objCnt + 1][d + dim] = cu[d];
		}

		Hypercube hc(dim, cl, cu);
		p[objCnt++] = new RtreeNodeEntry(id, hc);

		//log information
		if (objCnt % 1000 == 0)
			cout << ".";
		if (objCnt % 10000 == 0)
			cout << objCnt << " objects loaded" << endl;
	}

	double rawSize = dataSize(objCnt, dim);
	cout << "Total number of objects: " << objCnt << endl;
	_totalSpaceCost += rawSize;

	// build rtree
	cout << "Bulkloading R-tree..." << endl;
	const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
	mem = new FileMemory(PAGESIZE, indexfile, RtreeNodeEntry::fromMem, true);
	rtree = TGS::bulkload(*mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p, objCnt, false);

	cout << "[Rtree build done]" << endl;

	// load r-tree to memory, in-memory rtree
	cout << "loading R-tree for In-Memeory Processing..." << endl;
	//rtreeRAM(*rtree, ramTree);   
	_totalSpaceCost += ramTree.size() * 4096.00 / MB;
    cout << "Memory Size:" << totalSpaceCost << " MB" << endl;
	rtreeRAM(*rtree, ramTree);
}



int main(int argc,char **argv){
    clock_t start,finish;

    //Initialize 
    const char *datafile = param(argc,argv,"-f","data.txt");
    const char *indexfile = param(argc,argv,"-i","index.txt");
    const char *userfile = param(argc,argv,"-u","user.txt");
    int n = std::stoi(param(argc,argv,"-n","2")); // number of user
    int dim = std::stoi(param(argc,argv,"-d","2"));
    int k = std::stod(param(argc,argv,"-k","3"));
	m = std::stod(param(argc,argv,"-m","1"));
    int udis = std::stod(param(argc,argv,"-ud","1"));
    int pdis = std::stod(param(argc,argv,"-pd","1"));
    const char* method_name = param(argc, argv, "-method", "AA");
    load_user(user,dim,userfile,n);

    string ud = udis == 1 ? "IND" : "CL";
    string pd = "IND";
    if (pdis == 2){
	pd = "COR";
    }
    else if (pdis == 3){
	pd = "ANTI";
    }

    build_rtree(datafile,indexfile,dim,mem);

    string file_name = "./res/m_" +  to_string(m) + "_d_" + to_string(dim) + "_p_" + to_string(objCnt)
					+ pd + "_k_" + to_string(k) + "_u_" + to_string(user.size()) + ud + ".txt";
    string err_name = "./err/m_" +  to_string(m) + "_d_" + to_string(dim) + "_" + to_string(objCnt)
					 + "_k_" + to_string(k) + "_u_" + to_string(user.size()) + ".txt";

    freopen(file_name.c_str(),"w",stdout);
    freopen(err_name.c_str(),"w",stderr);

    printf("K: %d\n",k);
    printf("Data File : %s\n",datafile);
    printf("Index File: %s\n",indexfile);
    printf("Data Dimension : %d\n",dim);

	
    printf("User size: %d\n",user.size());
	

    cout << "Start solving...." << endl;

    start=clock();

    if (strcmp(method_name, "AA")) {
		Advanced solver(dim);
		solver.solve(k,m);	
    } else if (strcmp(method_name, "BSL")) {
		BaseLine solver(dim);
		sovler.solve(k,m);
    }
	

	
    finish=clock();

    printf("End solving\n");

    double totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"Total time : " << totaltime << "s" <<endl;
	
    int hyper_sum = 0;
    for (int i = 0 ; i < solver.m_final.size(); i++){
		hyper_sum += solver.m_final[i]->IntersectHP.size() +
					 solver.m_final[i]->AveHP.size() +
					 solver.m_final[i]->BelowHP.size();
    }

    cout << "Hyperplane size : " << hyper_sum << endl;
    cout << "Number of regions : " << num_cell << endl;
    return 0;
}
