#include "baseline.h"

using namespace std;
bool dcmp(double a,double b){
    return a > b;
}

struct Compare{
    bool operator()(const pair<double,Product>& x, const pair<double,Product>& y){
        return x.first > y.first ;
    }
};


BaseLine::BaseLine(const std::vector<User>& user,int dim)
    :m_user{user},
    m_dim{dim},
    m_celltree{nullptr}{

            
    }

BaseLine::~BaseLine(){
    delete m_celltree;
}

Product BaseLine::solve(int k,double ratio){
    gettopk(k);
    insert(ratio);
}

void BaseLine::gettopk(int k){

    queue<long> H;
    RtreeNode* node;

    Rtree rt = *rtree;
    
    long pageID;
    for (int i = 0 ; i < m_user.size(); i++){
        
        H.push(rt.m_memory.m_rootPageID);
        vector<double> scores;
        priority_queue<pair<double,Product>,vector<pair<double,Product>>,Compare> pq;

        while (!H.empty())
        {
            pageID = H.front();
            H.pop();
            node = rt.m_memory.loadPage(pageID);
            
            if (node->isLeaf() == false)
            {

                for (int j = 0; j < node->m_usedspace; j++)
                {   
                    if(pq.size() == k){
                        Hypercube& cube = node->m_entry[j]->m_hc;
                        const Point& center = cube.getUpper();
                        Product product(m_dim);
                        
             
                    
                       double score = m_user[i].score(center);
                       if(score < pq.top().first){
                           continue;
                       }
                    }
                
                    H.push(node->m_entry[j]->m_id);
                }
            }
            else {
                

                for (int j = 0 ; j < node->m_usedspace; j++){
                    

                    Hypercube& cube = node->m_entry[j]->m_hc; 
                    const Point& center = cube.getCenter();                   
                    Product product(m_dim);

                    for (int d = 0 ; d < m_dim; d++){         
                        product.set(d,center.m_coor[d]);
                    }
 
                    double score = m_user[i].score(product);
                    scores.push_back(score);
                    
                    if (pq.size() < k){
                        pq.push(make_pair(score,product));
                    }
                    else if(pq.size() == k){
                        if (score > pq.top().first){
                            pq.pop();
                            pq.push(make_pair(score,product));
                        }
                    }
                    
                
                }
                
            }

            
        } 
       // cout << "Score size : " << scores.size() << endl;
        sort(scores.begin(),scores.end(),dcmp);
     //   cout << "Score for user " << i << ":" << scores[k-1] << " " << pq.top().first << " " << pq.top().second << endl; 
        m_result.push_back(pq.top());
        //build halfspace for celltree
        vector<float> halfspace;
        for(int j = 0 ; j < m_dim; j++){
            halfspace.push_back(m_user[i].get(j));
        }
        halfspace.push_back((float)pq.top().first);
        HalfSpaces.push_back(halfspace);
        RecordIDtoHalfPlaneID[hsid.size()] = hsid.size();
        hsid.push_back(hsid.size()); //for plane id
        
    }
    
}

void BaseLine::insert(int m){
    if(m_celltree != nullptr){
        delete m_celltree;
    }
    Point pt;
    m_celltree = new cellTree();
    printf("Start insertion\n");
    m_celltree->insert(hsid,m,*rtree,pt,m_final);
    printf("End insertion\n");
    m_celltree->collectLeaf(m_final,m);
}

void BaseLine::optimize(){

}