#include "advanced.h"
#include <random>
#define DEBUG false

struct Triple{
    double score;
    long id;
    int pos;
    Triple(double _score,long _id,int _pos):
        score{_score},
        id{_id},
        pos{_pos}{
           
        }
};


struct Compare{
    bool operator()(const Triple& x, const Triple& y){
        return x.score > y.score ;
    }
};




Advanced::Advanced(int dim):Algorithm(),m_dim{dim},convex_sum{0}
        ,cover_sum{0},
        mbr_sum{0},
        insert_sum{0},
        pre_insert_sum{0},
        update_times{0},
        update_work{0}
        {

}
/**
* Calculate the k-th product of a user.
* @param index the index of user in global variable users
* @param k a positive integer 
* @return the index of product in R-tree
*/
pair<long,int> Advanced::topkitems(int index,int k){
    
    const Rtree& rt= *rtree;
    queue<long> H;
    long pageID;
    priority_queue<Triple,vector<Triple>,Compare> pq;
    RtreeNode* node;
    vector<float>& u = user[index];

    H.push(rt.m_memory.m_rootPageID);
    
    while(!H.empty()){
        pageID = H.front();
       
        H.pop();

        node = ramTree[pageID];

        if(node->isLeaf() == false){
         
            for (int j = 0; j < node->m_usedspace; j++)
                {   
                    if(pq.size() == k){
                       Hypercube& cube = node->m_entry[j]->m_hc;
                       const Point& upper = cube.getUpper();
                        
                       double score = user_score(u,upper);

                       if(score < pq.top().score){
                           continue;
                       }
                    }
                
                    H.push(node->m_entry[j]->m_id);
                }
        }
        else{
           
            for (int j = 0 ; j < node->m_usedspace; j++){
                    
                    
                    Hypercube& cube = node->m_entry[j]->m_hc; 
                    const Point& center = cube.getCenter();                   

 
                    double score = user_score(u,center);
                    
                    if (pq.size() < k){
                        pq.push(Triple(score,pageID,j));
                    }
                    else if(pq.size() == k){

                        if (score > pq.top().score){
                            pq.pop();
                            pq.push(Triple(score,pageID,j));
                        }
                    }
            }       
        }

       // avoid memory leak

    }
    pair<long,int> kth = make_pair(pq.top().id,pq.top().pos);
    
    vector<float> new_user = user[index];
    new_user.push_back(pq.top().score);
    HalfSpaces.push_back(new_user);
    hsid.push_back(index);

    return kth;
}

Advanced::~Advanced(){
    delete m_celltree;
}

/**
* The entry of advanced solution
* The result will be stored in m_final
* @param k a positive integer
* @param m a positive integer
* @return nothing
*/
void Advanced::solve(int k,int m){
    clock_t start_group = clock();
    this->grouping(k);
    clock_t end_group = clock();


    clock_t start_update ;
    clock_t end_update;
    double update_sum = 0;

    clock_t sum_start;
    clock_t sum_end;
    double sum_sum = 0;
    
    g_cnt = 0;
    g_prune = 0;

    int loop_times = 0;
    if (m_celltree != nullptr){
        delete m_celltree;
    }

    m_celltree = new cellTree();
    
    for (int i = 0 ; i < m_group.size(); i++){
        m_celltree->root->G.push_back(i);
    }
    pq.push(m_celltree->root);
    // cell's bound
    for (int i = 0; i < m_dim; i++){

        for (int j = 0 ; j < 2; j++){
            vector<float> bound(m_dim+1,0);
            bound[i] = 1;
            bound[m_dim] = j;
            HalfSpaces.push_back(bound);
            if (j == 0){
                m_celltree->root->IntersectHP[HalfSpaces.size() - 1] = true;
            }
            else{
               m_celltree->root->IntersectHP[HalfSpaces.size() - 1] = false;
            }
        }
        
    }



    m_celltree->root->upper = user.size();
    m_celltree->root->depth = 0;
    bool is_root = true;
   
    while (!pq.empty()){
       
        loop_times++;
   
        cell* c = pq.top();
       
        pq.pop();   
      
        if (is_root == false){
           start_update = clock();
          
           int is_optimal;
           is_optimal = update_group_bo(c,m);
         
            
           end_update = clock();
           update_sum += (end_update - start_update) / CLOCKS_PER_SEC;
          
           if (is_optimal == 2){
                // if lp_solve returns 2, it means the cell is infesiable
               continue;
           }
            
           if (c->rank >= m){ //c->lower + sum >= m cannot be used for pruning whole cell
                m_final.push_back(c);
                continue;
            }

            if ( c->upper < m || c->G.size()==0){
                continue;
            }
        }

        int index = 0;
       
         // get minimal group
        if (is_root == false){
           
            //int minval = 0;
            int minval = 0x7fffffff;
            for (int i = 0; i < c->G.size(); i++){
                //if (m_group[c->G[i]].user.size() > minval){
                if (m_group[c->G[i]].user.size() < minval){
                    minval = m_group[c->G[i]].user.size();
                    index = i;
                }
            }
        }
       
        else{
            int minval = 100;  // identify closest group

            for (int i = 0; i < c->G.size(); i++){
                RtreeNode* entry = ramTree[m_group[c->G[i]].kth.first];
                const Point& p = entry->m_entry[m_group[c->G[i]].kth.second]->m_hc.getCenter();
                double dis = 0;
                for (int i = 0; i < p.m_dimen; i++){
                    dis += (1 - p[i]) * (1 - p[i]);
                }
                if (dis < minval){
                    minval = dis;
                    index = i;
                }
            }
        }
          clock_t insert_start = clock();
          vector<cell*> leaves = insert(c,index);
           clock_t insert_end = clock();

         insert_sum += 1.0 * (insert_end - insert_start) / CLOCKS_PER_SEC;
        for (int i = 0; i < leaves.size(); i++){
         
            if (leaves[i]->rank >= m){
                m_final.push_back(leaves[i]);
            }
            else if(leaves[i]->upper <m){
                continue;
            }
            else if(leaves[i]->upper >= m && leaves[i]->G.size() > 0){
                pq.push(leaves[i]);
            }
        }
        is_root = false;

    }
    cout << "Update time : " << update_times << " Work : " << update_work << endl;
    cout << "Grouping time : " << (end_group - start_group) / CLOCKS_PER_SEC << endl;
    cout << "Inserting time : " << insert_sum << endl ;
   
    cout << "Convex time :" << convex_sum << endl;
    cout << "Cover time : " << cover_sum << endl;
    cout << "MBR time : " << mbr_sum << endl;
    cout << "Update time : " << update_sum << endl;
    cout << "Loop times : " << loop_times << endl;
    cout << "Number of leaves : " << m_final.size() << endl;

} 
/**
* Calculate the convex hull of points
* @param points a vector store all points 
* @return the indices of points which construct the convex hull
*/
vector<int> Advanced::convex_hull(vector<vector<float>>& points){
    using namespace orgQhull;

    using namespace std;

  
    clock_t convex_start = clock();
    vector<int> result;
    
   
    if (m_dim >= 3){
        RboxPoints rbox;
        Qhull qhull;
        string s = to_string(m_dim - 1) + " " + to_string(points.size()) + " ";
        for (int i = 0; i < points.size(); i++){
            for (int j = 0; j < m_dim - 1; j++){
                s += to_string(points[i][j]) + " ";
            }
            
        }
   

        std::istringstream is(s);

        rbox.appendPoints(is);

        qhull.runQhull(rbox,"");

        std::stringstream output;
        qhull.setOutputStream(&output);
        qhull.outputQhull("Fx");

        int num = 0;
        output >> num;
       
        int tmp;
        for (int i = 0; i < min(num,(int)points.size()); i++){
            output >> tmp;
            result.push_back(tmp);
        }

     }
     else { // for 2-d case
        int min_w = 0x7fffffff;
        int min_w_idx = 0;
        int max_w = -1;
        int max_w_idx = 0;
        for (int i = 0; i < points.size(); i++){
            if (points[i][0] > max_w){
                max_w = points[i][0];
                max_w_idx = i;
            }

            if (points[i][0] < min_w){
                min_w = points[i][0];
                min_w_idx = i;
            }
        }

        result.push_back(min_w_idx);
        result.push_back(max_w_idx);
     }
    clock_t convex_end = clock();
    convex_sum += 1.0 * (convex_end - convex_start) / CLOCKS_PER_SEC;
    return result;
}
/**
* Insert the hyperplane
* @param c a pointer points to a cell
* @param index the index of group to be inserted in cell c
* @return new cells generated by the insertion of hyperplanes
*/
vector<cell*> Advanced::insert(cell *c, int index){
  
 
    return insert_chull(c,index); 
   
}


/**
* Insert the convex hull of a user group in a cell 
* @param c a pointer points to a cell
* @param index the index of group to be inserted in cell c
* @return new cells generated by the insertion of hyperplanes
*/
vector<cell*> Advanced::insert_chull(cell* c, int index){
    
    vector<cell*> result;
    
    if (c == m_celltree->root){
        clock_t start_mbr = clock();
        lprec* lp = make_lp(0, m_dim);
        m_celltree->lpModel(lp, m_dim, c->IntersectHP);
        vector<float> cl, cu;
        int is_optimal =  m_celltree->findCellMBR(lp, cl, cu);
        c->cl = cl;
        c->cu = cu;
        delete_lp(lp);
    }
    vector<vector<float>> vertices = get_vertices(c->cl,c->cu);
    

    vector<long> insert_hyper;
    // if group size is large enough for finding a convex hull
    if (m_group[c->G[index]].user.size() >= m_dim ){

        vector<vector<float>> points;
        for (int i = 0 ; i < m_group[c->G[index]].user.size() ; i++){
            points.push_back(user[m_group[c->G[index]].user[i]]);
        }
        vector<int> chull_idx = convex_hull(points);    
        vector<int> cover, inter, uncover ,hypers;

        for (int i = 0; i < chull_idx.size(); i++){
            hypers.push_back(m_group[c->G[index]].user[chull_idx[i]]);
        }

        cover_check(c, hypers, vertices, cover, inter, uncover);

        c->rank += cover.size();
        c->lower += cover.size();
        c->upper -= uncover.size();
        if (c->rank >= m){
            m_final.push_back(c);
            return result;
        }
        else if (c->upper < m){
            return result;
        }

        int insert_limit = 8;
        if (c->depth > 10){
            insert_limit = 5;
        }
        else if (c->depth > 20){
            insert_limit = 3;
        }
        cout << inter.size() << endl;
        for (int i = 0; i < inter.size() && i < insert_limit; i++){
          
            
            insert_hyper.push_back(m_group[c->G[index]].user[chull_idx[inter[i]]]);
            
        }
      
        Group backup = m_group[c->G[index]];
        c->G.erase(c->G.begin() + index);
        
        for(int i=0; i<cover.size(); i++)
        {
            vector<int>::iterator iter = find(backup.user.begin(),backup.user.end(),cover[i]);
            if(iter!=backup.user.end())
            {
                backup.user.erase(iter);
            }
        }
        for(int i=0; i<uncover.size(); i++)
        {
            vector<int>::iterator iter = find(backup.user.begin(),backup.user.end(),uncover[i]);
            if(iter!=backup.user.end())
            {
                backup.user.erase(iter);
            }
        }
        for(int i=0; i<insert_hyper.size(); i++)
        {
            vector<int>::iterator iter = find(backup.user.begin(),backup.user.end(),insert_hyper[i]);
            if(iter!=backup.user.end())
            {
                backup.user.erase(iter);
            }
        }

        if (backup.user.size() > 0){
            c->G.push_back(m_group.size());
            m_group.push_back(backup);
        }
    }
    else{ // group size is small
        vector<int> cover, inter, uncover;
        cover_check(c,m_group[c->G[index]].user,vertices, cover, inter, uncover);
        c->rank += cover.size();
        c->lower += cover.size();
        if (c->rank >= m){
            m_final.push_back(c);
            return result;
        }
        c->upper -= uncover.size();
        if (c->upper < m){
            return result;
        }
        
        for (int i = 0; i < inter.size(); i++){
            insert_hyper.push_back(m_group[c->G[index]].user[inter[i]]);
        }

        c->G.erase(c->G.begin() + index);
        
        
    }

   int count = 0;
    
    if (insert_hyper.size() > 0){
        m_celltree->insert(c, insert_hyper, m, result);
    }
    return result;
}


/**
* Compare the pivot in the insertion
* @param lower the lower bound point of a cell
* @param higher a product point
* @return true if the lower bound dominates the product point
*		  false if the  product point dominates the lower bound 
*/
bool Advanced::comparepivot(vector<float>& lower, const Point& higher)
{
    for(int i=0; i<lower.size(); i++)
    {
        if(lower[i] > higher.m_coor[i])
            return false;
    }
    return true;
}

/**
* Compare the pivot in the insertion
* @param lower the upper bound point of a cell
* @param a product point
* @return true if the
*/
bool Advanced::comparepivotII(const Point& lower, vector<float>& higher)
{
    for(int i=0; i<higher.size(); i++)
    {
        if( lower.m_coor[i] > higher[i])
            return false;
    }
    return true;
}

/**
* Update the group before insertion
* @param c a pointer points to cell 
* @param m a positive integer 
* @return 0 if the update is successful
*         1 if the update is falied (mbr calculation)
*/
int Advanced::update_group_bo(cell *c,int m){
  
    // get MBR (lower and upper) of current cell
    clock_t start_mbr = clock();
    lprec* lp = make_lp(0, m_dim);
    m_celltree->lpModel(lp, m_dim, c->IntersectHP);
    vector<float> cl, cu;
    int is_optimal =  m_celltree->findCellMBR(lp, cl, cu);
    delete_lp(lp);
    // get vertices of mbr


    clock_t end_mbr = clock();
    mbr_sum += 1.0 * (end_mbr - start_mbr) / CLOCKS_PER_SEC;

    if (is_optimal == 2){
        return is_optimal; // if infeasible , ignore this cell
    }

    c->cu = cu;
    c->cl = cl;

    vector<vector<float>> vertices = get_vertices(cl,cu);
    vector<cell*> cv;
    unordered_map<int,vector<int>> cover, uncover;
    vector<int> removegroups;

    // check groups first
    for (vector<int>::iterator iter = c->G.begin(); iter != c->G.end();){
        update_times++;
        RtreeNode* node = ramTree[m_group[*iter].kth.first];
        Hypercube& cube = node->m_entry[m_group[*iter].kth.second]->m_hc; 
        const Point& kthpoint = cube.getCenter();     
        if(comparepivot(cl,kthpoint))
        {
            update_work++;
            c->rank += m_group[*iter].user.size();
            if(c->rank >= m)
            {
                m_final.push_back(c);
                return 0;
            }
            c->G.erase(iter);
            continue;
        }
        else if(comparepivotII(kthpoint,cu))
        {
            update_work++;
            c->upper -= m_group[*iter].user.size();
            if(c->upper < m)
            {
                return 0;
            }
            c->G.erase(iter);
            continue;
        }
        iter++;
    }

    for(vector<int>::iterator iter = c->G.begin(); iter != c->G.end(); iter++){
         // get covnvex hull
        update_times++;
        vector<vector<float>> users;
        // index to user
   
        // get hyperplanes of users ( convex hull vertices)
        vector<int> hypers;
        if (m_group[*iter].user.size() >= m_dim){
             for (int i = 0; i < m_group[*iter].user.size(); i++){
                users.push_back(user[m_group[*iter].user[i]]);
            }
            vector<int> cH = convex_hull(users);
            for (int i = 0; i < cH.size(); i++){
                hypers.push_back(m_group[*iter].user[cH[i]]);
            }
        }
        else {
            hypers = m_group[*iter].user;
        }
        //group_cover_check(c, *iter, hypers, removegroups, vertices, cover, uncover);       
        bool flag = true; // full or partial indicator 
        Group groupCopy = m_group[*iter];
        group_cover_check_full(flag, c, *iter, groupCopy, hypers, removegroups, vertices, cover, uncover);
    }
    if(c->rank<m&&c->upper>m)
    {
        erase_pruned_groups_users(c, removegroups, cover, uncover);
    }
    return 0;
        
        
}



/**
* Group all users based on their k-th product
* Add group information into m_group
* @param k a positive integer
* @return nothing
*/
void Advanced::grouping(int k){

    for (int i = 0; i < user.size(); i++){
        if (i % 1000 == 0)
            cout << i << endl;
        pair<long,int> res = this->topkitems(i,k);
        
        bool flag = false;
        
        for (int j = 0 ;j < m_group.size(); j++){
            Group& g = m_group[j];
               
            if (res == g.kth){  
                g.user.push_back(i);
                flag = true;
                break;
            }
        }

        if (flag == false){
           
           Group g;
      
           g.kth = res;
           g.user.push_back(i);
           m_group.push_back(g);

        }
    }
   
    
   
    int max_val = -1;
    int min_val = 100000;
    for (int i = 0 ;i < m_group.size(); i++){
        max_val = max(max_val,(int)m_group[i].user.size());
        min_val = min(min_val,(int)m_group[i].user.size());

     
        RtreeNode* node = rtree->m_memory.loadPage(m_group[i].kth.first);
        const Hypercube& cube = node->m_entry[m_group[i].kth.second]->m_hc;
        const Point& p = cube.getCenter();
     

        delete node;
        
    } 
    cout << "Number of group : " << m_group.size() << endl;
    cout << "Max group size :" << max_val << endl;
    cout << "Min group size :" << min_val << endl;
  

    
}

/**
* Generate the hyperplane of user
* @param option user's k-th product
* @param user user preference
* @return the user's hyperplane in product space
*/
vector<float> Advanced::hyperplane(const Point& option, vector<float>& user){
    double score = 0.0;
    vector<float> res(user);
    
    for (int i = 0 ; i < user.size() ; i++){
        score += user[i] * option.m_coor[i];
    }
    res.push_back(score);
    return res;
}

/**
* Generate the 0-1 string for MBR points
* @param depth a positive integer  the length of string
* @param gen generated strings
* @param seq 0-1 integer sequence in dfs
* @return nothing
*/
void Advanced::gen_string(int depth,vector<string>& gen,vector<int>& seq){
    if (depth == m_dim){
        string s = "";
        for(int i = 0 ; i < seq.size(); i++){
            if (seq[i] == 1){
                s += "1";
            }
            else{
                s += "0";
            }
           
        }
        
        gen.push_back(s);
     
        return ;
    }
      
    for (int i = 0 ; i < 2; i++){
       
        seq.push_back(i);
        gen_string(depth+1,gen,seq);
        seq.pop_back();
        
    }
}
/**
* Generate vertices of MBR 
* @param cl lower bound point of MBR
* @param cu upper bound point of MBR
* @return 2^d points (d is the dimension of upper bound point and lower bound point)
*/
vector<vector<float>> Advanced::get_vertices(vector<float>& cl, vector<float>& cu){
   
    vector<vector<float>> result;
    vector<string> strs;
    vector<int> seq;
    gen_string(0,strs,seq);
   
    for (int i = 0; i < strs.size(); i++){
        vector<float> item;
      
        for (int j =0; j < strs[i].size(); j++){
            if (strs[i][j] == '1'){
                item.push_back(cu[j]);
            }
            else{
                item.push_back(cl[j]);
            }
        }
        result.push_back(item);   
    }
    

    return result;
}
/**
* Check the realtionship between the cell and hyperplanes
* @param c a pointer points to a cell
* @param hypers indices of hyperplanes in global variable HalfSpaces
* @param vertices vertices of MBR
* @param cover indices of hyperplane which cover the cell
* @param inter indices of hyperplane which intersects the cell
* @param uncover indices of hyperplane which uncover the cell
* @return nothing
*/
void Advanced::cover_check(cell *c,vector<int>& hypers, vector<vector<float>>& vertices, vector<int>& cover, vector<int>& inter, vector<int>& uncover){
    unordered_map<long int, bool> touchhs;
    for (auto iter = c->IntersectHP.begin(); iter != c->IntersectHP.end(); iter++)
    {
        touchhs[iter->first] = iter->second;
    }  
    inter.clear();
    cover.clear();
    uncover.clear();
    clock_t cover_start = clock();
    for (int i = 0; i < hypers.size(); i++){
       
        int cnt = 0;
        int ncnt = 0;
        vector<float>& hyper = HalfSpaces[hypers[i]];   
        for (int j = 0 ;j < vertices.size(); j++){
            float sum = 0;
            for (int k = 0 ; k < vertices[j].size(); k++){
                sum += hyper[k] * vertices[j][k];
            }
            
            if ( sum >= hyper[hyper.size() - 1]){
               cnt++;
            }
            else
            {
                ncnt++;
            }
            if(cnt>0 && ncnt >0)
            {
                break;
            }
        }

        if (cnt == vertices.size()){
            cover.push_back(i);
         
        }
        else if (ncnt == vertices.size()){
            uncover.push_back(i);          
        }
        else 
        {
            if (m_celltree->isFeasible(touchhs,hypers[i],false) == false){
            //cover 
                cover.push_back(i);
            }
            else if ( m_celltree->isFeasible(touchhs,hypers[i],true) == false){
            //uncover
                   uncover.push_back(i);
            }
            else{
                    inter.push_back(i);
            }
        }
    }

    if (cover.size()  + uncover.size() + inter.size() != hypers.size()){
        cout << "Check error " << cover.size()  << " " << uncover.size() << " " << inter.size() << " "<<hypers.size() << endl;
        exit(-1);
    }
    clock_t cover_end = clock();

    cover_sum += 1.0 * ( cover_end - cover_start ) / CLOCKS_PER_SEC; 
        
}


/**
* Erase pruned groups and users
* @param c a pointer points to a cell
* @param removegroups groups to be pruned
* @param cover users cover the cell
* @param uncover users uncover the cell
* @return nothing
*/
void Advanced::erase_pruned_groups_users(cell* c, vector<int>& removegroups, unordered_map<int,vector<int>>& covers, unordered_map<int, vector<int>>& uncover)
{


    unordered_map<int, Group> updateuGroup;
    for(unordered_map<int, vector<int>>::iterator iter = covers.begin(); iter!=covers.end(); iter++)
    {
        if(updateuGroup.find(iter->first)==updateuGroup.end())
        {
            updateuGroup[iter->first] = m_group[iter->first];
        }
        for(int ruid=0; ruid<iter->second.size(); ruid++)
        {
            vector<int>::iterator uiter = find(updateuGroup[iter->first].user.begin(), updateuGroup[iter->first].user.end(), iter->second[ruid]);
            if(uiter !=updateuGroup[iter->first].user.end())
            {
                updateuGroup[iter->first].user.erase(uiter);
            }
        }
    }

    // remove items in uncover group
    for(unordered_map<int, vector<int>>::iterator iter = uncover.begin(); iter!=uncover.end(); iter++)
    {
        if(updateuGroup.find(iter->first)==updateuGroup.end())
        {
            updateuGroup[iter->first] = m_group[iter->first];
        }
        for(int ruid=0; ruid<iter->second.size(); ruid++)
        {
            vector<int>::iterator uiter = find(updateuGroup[iter->first].user.begin(), updateuGroup[iter->first].user.end(), iter->second[ruid]);
            if(uiter !=updateuGroup[iter->first].user.end())
            {
                updateuGroup[iter->first].user.erase(uiter);
            }
        }
    }

    // remove these groups which can be pruned definitly.
    
    for(int rgid=0; rgid <removegroups.size(); rgid++)
    {
        update_work += removegroups.size();
        vector<int>::iterator iter = find(c->G.begin(), c->G.end(), removegroups[rgid]);
        if (iter != c->G.end())
            c->G.erase(iter);
    }

    // replace group by the updated partial group
    for(unordered_map<int, Group>::iterator giter = updateuGroup.begin(); giter!=updateuGroup.end(); giter++)
    {
        vector<int>::iterator giditer = find(c->G.begin(), c->G.end(), giter->first);
        if (giditer != c->G.end())
            c->G.erase(giditer);
        c->G.push_back(m_group.size());
        m_group.push_back(giter->second);
    }
}

/**
* Group check in updating
* @param flag a bool flag
* @param c a pointer points to cell
* @param id group to be checked
* @param group_copy checked group
* @param hypyers user hyperplanes to be checked
* @param removegroups groups to be pruned (result)
* @param vertices vertices of MBR
* @param cover users cover the cell (result)
* @param uncover users uncover the cell (result)
* @return nothing
*/
void Advanced::group_cover_check_full(bool& flag, cell *c, int& id, Group& group_copy, 
vector<int>& hypers, vector<int>& removegroups, vector<vector<float>>& vertices, 
unordered_map<int,vector<int>>& cover, unordered_map<int,vector<int>>& uncover)
{
    clock_t cover_start = clock();

    unordered_map<long int, bool> touchhs;
    for (auto iter = c->IntersectHP.begin(); iter != c->IntersectHP.end(); iter++)
    {
        touchhs[iter->first] = iter->second;
    }  

    // full cover
    unordered_map<long int, bool> fullcoverhs = touchhs;
    for (int i = 0; i < hypers.size(); i++)
    {
        fullcoverhs[hypers[i]] = true;
    }
    if(m_celltree->isCellFeasible(fullcoverhs)==true)
    {
        c->rank += group_copy.user.size();
        if(c->rank >= m)
        {
            m_final.push_back(c);
            return;
        }
        for(int hid=0; hid <group_copy.user.size(); hid++)
        {
            c->IntersectHP[group_copy.user[hid]] = true;
        }
        if(flag==true) // the whole group
            removegroups.push_back(id);
        else // partial of the group
        {
            if(cover.find(id)==cover.end())
            {
                cover[id] = group_copy.user;
            }
            else
            {
                cover[id].insert(cover[id].end(), group_copy.user.begin(), group_copy.user.end());
            }
        }
    }
    // non cover
    unordered_map<long int, bool> uncoverhs = touchhs;
    for (int i = 0; i < hypers.size(); i++)
    {
        uncoverhs[hypers[i]] = false;
    }
    if(m_celltree->isCellFeasible(uncoverhs)==true)
    {
        c->upper -= group_copy.user.size();
        if(c->upper < m)
        {
            return;
        }
        if(flag==true) // the whole group
            removegroups.push_back(id);
        else // partial of the group
        {
            if(uncover.find(id)==uncover.end())
            {
                uncover[id] = group_copy.user;
            }
            else
            {
                uncover[id].insert(uncover[id].end(), group_copy.user.begin(), group_copy.user.end());
            }
            
        }
        
    }

    
    // intersects
    for (int i = 0; i < hypers.size(); i++)
    {   
        int cnt = 0;
        int ncnt = 0;
        vector<float>& hyper = HalfSpaces[hypers[i]];   
        for (int j = 0 ;j < vertices.size(); j++){
            float sum = 0;
            for (int k = 0 ; k < vertices[j].size(); k++){
                sum += hyper[k] * vertices[j][k];
            }
            
            if ( sum >= hyper[hyper.size() - 1]){
               cnt++;
            }
            else
            {
                ncnt++;
            }
            if(cnt>0 && ncnt >0)
            {
                break;
            }
        }

        if (cnt == vertices.size()){
            cover[id].push_back(i);
            c->rank++;
            if(c->rank >= m)
            {
                m_final.push_back(c);
                return;
            }
        }
        else if (ncnt == vertices.size()){
            uncover[id].push_back(i);          
            c->upper--;
            if(c->upper < m)
            {
                return;
            }
        }
        else 
        {
            if (m_celltree->isFeasible(touchhs,hypers[i],false) == false){
            //cover 
                cover[id].push_back(i);
                c->rank++;
                if(c->rank >= m)
                {
                    m_final.push_back(c);
                    return;
                }
            }
            else if ( m_celltree->isFeasible(touchhs,hypers[i],true) == false){
            //uncover
                uncover[id].push_back(i);
                c->upper--;
                if(c->upper < m)
                {
                    return;
                }
            }
        }
    }

    // remove convex hull users from goupidcopy 
    for(int userid=0; userid < hypers.size(); userid++)
    {
        vector<int>::iterator iter = find(group_copy.user.begin(), group_copy.user.end(), hypers[userid]);
        if (iter != group_copy.user.end())
            group_copy.user.erase(iter);
    }
    hypers.clear();
    
   
    clock_t cover_end = clock();

    cover_sum += 1.0 * ( cover_end - cover_start ) / CLOCKS_PER_SEC; 
        
}

/**
* Check the realtionship between the cell and hyperplanes without MBR check 
* @param c a pointer points to a cell
* @param hypers indices of hyperplanes in global variable HalfSpaces
* @param vertices vertices of MBR
* @param cover indices of hyperplane which cover the cell
* @param inter indices of hyperplane which intersects the cell
* @param uncover indices of hyperplane which uncover the cell
* @return nothing
*/
void Advanced::exactly_check(vector<int>& hypers,cell* c,vector<int>& cover, vector<int>& inter, vector<int>& uncover){
    
    clock_t cover_start = clock();
   
    unordered_map<long int, bool> touchhs;
    for (auto iter = c->IntersectHP.begin(); iter != c->IntersectHP.end(); iter++)
    {
        touchhs[iter->first] = iter->second;
    }  
    for (int i = 0; i < hypers.size(); i++){
        
        if (m_celltree->isFeasible(touchhs,hypers[i],false) == false){
            //cover 
            cover.push_back(i);
        }
        else if ( m_celltree->isFeasible(touchhs,hypers[i],true) == false){
            //uncover
            uncover.push_back(i);
        }
        else{
            inter.push_back(i);
        }

    } clock_t cover_end = clock();

    cover_sum += 1.0 * ( cover_end - cover_start ) / CLOCKS_PER_SEC;     
}