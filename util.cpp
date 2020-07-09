#include "util.h"

using namespace std;

char* param(int argc,char **argv,char *target,char *def){
    for(int i = 1 ; i < argc - 1; i++){
        if(strcmp(argv[i],target) == 0){
            return  argv[i+1];
        }
    }
    return def;
}


void load_user(vector<vector<float>>& user,int dim,const char* path,int n){
    
    fstream file(path);
    
    double d;
    int temp;
    for (int i = 0 ; i < n ;i++){

        file >> temp ;
        vector<float> u;
        float rect[dim * 2 + 1];
        for (int i = 0 ; i < dim * 2; i++){
            file >> rect[i];

        } 
     
        for (int i = 0 ; i < dim ; i++){
            u.push_back( (rect[i] + rect[i + dim] ) / 2.0);
        }
         user.push_back(u);
      
       
    }

    
}
bool cmp(double a,double b){
    return a > b;
}

double user_score(const vector<float>& user,const Point& p){
    
    double score = 0.0;

    for (int i = 0 ; i < user.size(); i++){
        
        score += user[i] * p.m_coor[i];
        
    }
    
    return score;
}

double user_score(const vector<float>& user,const vector<float>& product){
    
    double score = 0.0;

    for (int i = 0 ; i < user.size(); i++){
        
        score += user[i] * product[i];
        
    }
    
    return score;
}


void testtopk(const char* datafile,int k,int dim,int psize,const vector<vector<float>>& user){
    fstream df(datafile);

    int temp;
    vector<vector<float>> product;

    for (int i = 0; i < psize ; i++){
        double l[dim];
        double u[dim];

        df >> temp; 
     
       for (int j = 0 ; j < dim; j++){
            df >> l[j];

        }
        for (int j = 0 ;j < dim ; j++){
            df >> u[j];
        }
        
        vector<float> p;
        for (int j = 0; j < dim ;j++){
            p.push_back((l[j] + u[j]) / 2.0);  

        }

        product.push_back(p);
       
    }
    
   
     for (int i = 0 ; i < user.size(); i++){

        vector<double> score;
        
        for (int j = 0 ;j < product.size(); j++){
            score.push_back(user_score(user[i],product[j]));
        }

        sort(score.begin(),score.end(),cmp);
        cout << "Score for user " << i << " :" << score[k-1] << endl;
    }
}

void visit_hs(int id){

	cout << "(";
	for(int i = 0 ; i < HalfSpaces[id].size(); i++){
		cout << HalfSpaces[id][i] << " ";
	}
	cout << ")" << endl;
}

void display_cell(std::vector<cell*>& cells){
    cout << "==============================" << endl;
	cout << "Leaves :" << cells.size() << endl;
	for(int i = 0 ; i < cells.size(); i++){

		cout << "--------------------" << endl;
		cout << "Cell " << i << endl;
		
		cell& c = *cells[i];
		cout << "count :" << c.rank << endl;
		cout << "above:" << endl;
		for(auto iter = c.AboveHP.begin(); iter != c.AboveHP.end() ; iter++){
			visit_hs(*iter);
		}
	
		cout << "below:" << endl;
		for(auto iter = c.BelowHP.begin(); iter != c.BelowHP.end() ; iter++){
			visit_hs(*iter);
		}
		cout << "intersect:" << endl;
		for(auto iter = c.IntersectHP.begin(); iter != c.IntersectHP.end() ; iter++){
			visit_hs( (*iter).first);
			cout << "Inter :" << (*iter).second << endl;
		}
		cout << "--------------------" << endl;

	}
}


void visit_rtree(Rtree& rt){
    ramTree.clear();
	queue<long> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long pageID;
	while (!H.empty())
	{
		pageID = H.front();
        cout << pageID << endl;
		H.pop();
		node = rt.m_memory.loadPage(pageID);
		
    	if (node->isLeaf() == false)
		{

            cout << "Not leaf" << endl;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
        else {
            
            cout << "Space :" << node->m_usedspace << endl;
			objCnt += node->m_usedspace;

            for (int i = 0 ; i < node->m_usedspace; i++){
				

                Hypercube& cube = node->m_entry[i]->m_hc;
                const Point& lower = cube.getLower();
                const Point& higher = cube.getUpper();
                cout << endl;
                cout << "Lower: " ;
                for(int i = 0 ; i < lower.m_dimen; i++){
                    cout << lower.m_coor[i] << " ";
                }
                cout << endl;

                cout << "Higher: " ;
                for(int i = 0 ; i < higher.m_dimen; i++){
                     cout << higher.m_coor[i] << " ";
                }
             
            }
            
        }

		
	}
}