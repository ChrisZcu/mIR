#include "cellTree.h"

int remain; //for pruning rule 2
int early_pruning;


bool dominateRecords(float* a, float* b, const int dimen)
{
	for (int i = 0; i < dimen; i++)
	{
		if (a[i] < b[i])
			return false;
	}
	return true;
}

cellTree::cellTree()
{
	root = new cell();
}

cellTree::~cellTree()
{
	for (auto iter = dagNode.begin(); iter != dagNode.end(); iter++)
	{
		iter->second->rDominator.clear();
		unordered_set<long int>().swap(iter->second->rDominator);
		delete iter->second;
	}
	dagNode.clear();
	unordered_map<long int, gNode*>().swap(dagNode);
	//releaseCell(root);
}

void cellTree::releaseCell(cell* node)
{
	cell* tmp;
	for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);

	while (node != nullptr)
	{
		for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);
		cell* old = node;
		node = node->left;
		old->IntersectHP.clear();
		unordered_map<long int, bool>().swap(old->IntersectHP);
		old->BelowHP.clear();
		unordered_set<long int>().swap(old->BelowHP);
		old->AboveHP.clear();
		unordered_set<long int>().swap(old->AboveHP);
		delete old;
	}
}

void cellTree::lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs)
{
	double row[MAXDIMEN];
	if (lp == NULL)
	{
		fprintf(stderr, "Unable to create new LP model\n");
		exit(0);
	}

	// constraints in each dimension
	for (int di = 0; di < dimen; di++)
	{
		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			if (dii - 1 == di)
			{
				row[dii] = 1;
			}
			else
			{
				row[dii] = 0;
			}
		}

		add_constraint(lp, row, GE, 0);
		add_constraint(lp, row, LE, 0x7fffffff);
	}

	// in reduced space, sum_{q_i} should less than 1
	/*for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	add_constraint(lp, row, GE, 0);
	add_constraint(lp, row, LE, 1);*/

	// constraints in intersected hyperplanes
	for (unordered_map<long int, bool>::iterator iter = touchhs.begin(); iter != touchhs.end(); iter++)
	{
		addHP(lp, iter->first, iter->second);
	}

	// set scale 
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void cellTree::addHP(lprec* model, long int hpid, bool sideindicator)
{
	int hsID = hpid;
	int dimen = HalfSpaces[0].size() - 1;
	double row[MAXDIMEN];

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = HalfSpaces[hsID][dii - 1];
	}
	if (sideindicator == false)
	{
		add_constraint(model, row, LE, HalfSpaces[hsID][dimen]);
	}
	else if (sideindicator == true)
	{
		add_constraint(model, row, GE, HalfSpaces[hsID][dimen]);
	}
	else
	{
		std::cout << "Unable to detect half plane direction!!!" << endl;
	}
}

bool cellTree::isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator)
{	
	num_lp++;
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 1;
	lprec *lp = make_lp(0, dimen);
	lpModel(lp, dimen, touchhs);
	addHP(lp, hpid, sideindicator);

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = -1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int ret = solve(lp);
	get_variables(lp, row);
	delete_lp(lp);

	if (ret == 0)
		return true;
	else
		return false;
	return true;
}

bool cellTree::isCellFeasible(unordered_map<long int, bool>& touchhs)
{	
	num_lp++;
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 1;
	lprec *lp = make_lp(0, dimen);
	lpModel(lp, dimen, touchhs);

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = -1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int ret = solve(lp);
	get_variables(lp, row);
	delete_lp(lp);

	if (ret == 0)
		return true;
	else
		return false;
	return true;
}

void cellTree::updateRank(cell* node, const int mink)
{	
	/*
	map<cell*,int>::iterator iter = this->inserted.find(node);

	if (iter == this->inserted.end()){
		this->inserted[node] = 1;
	}
	else{
		this->inserted[node]++;
	}
	if (inserted[node] == this->hyper_size){
		node->rank += this->chull_size - inserted[node];
		node->lower += this->chull_size - inserted[node];


	}
	else{
		
		
		
	}*/
	
	node->rank++;
	node->lower += 1;
	if (node->rank >= mink)
	{
		node->isPruned = true;
		if (node->left != NULL)
		{
			releaseCell(node->left);
			node->left = NULL;
		}
		if (node->right != NULL)
		{
			releaseCell(node->right);
			node->right = NULL;
		}

	}
	else
	{
		if (node->left != NULL)
			updateRank(node->left, mink);
		if (node->right != NULL)
			updateRank(node->right, mink);
	}
}

void cellTree::updateUpper(cell* node, const int mink)
{	
	node->upper--;
	if (node->upper < mink)
	{
		node->isPruned = true;
		if (node->left != NULL)
		{
			releaseCell(node->left);
			node->left = NULL;
		}
		if (node->right != NULL)
		{
			releaseCell(node->right);
			node->right = NULL;
		}

	}
	else
	{
		if (node->left != NULL)
			updateUpper(node->left, mink);
		if (node->right != NULL)
			updateUpper(node->right, mink);
	}
}

void cellTree::inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs)
{
	if (node->isPruned == false)
	{	
		for (auto iter = node->IntersectHP.begin(); iter != node->IntersectHP.end(); iter++)
		{
			touchhs[iter->first] = iter->second;
		}
		if (isFeasible(touchhs, hpid, false) == false)
		{
		//	node->BelowHP.insert(hpid);
			updateRank(node, mink);
			totalSpaceCost += 4;
		}
		else if (isFeasible(touchhs, hpid, true) == false)
		{
			updateUpper(node, mink);
		//	node->AboveHP.insert(hpid);
			totalSpaceCost += 4;
		}
		else
		{	
			
			if (node->left == NULL&&node->right == NULL)
			{
				
 				if (node->rank >= mink || node->upper < mink)
				{
					node->isPruned = true;
				}
			
				else
				{
					node->left = new cell(node);
					node->left->IntersectHP[hpid] = false; // above / uncover
					updateUpper(node->left, mink);
					if (node->left->upper < mink)
					{
						node->left->isPruned = true;
					}		

					node->right = new cell(node);
					node->right->IntersectHP[hpid] = true; // below / cover
 					updateRank(node->right, mink);

					if (node->right->rank >= mink)
					{
						node->right->isPruned = true;
					}		
				}
			}
			else if (node->left != NULL&&node->right != NULL)
			{	
				
				if (node->left->isPruned == false)
				{
					unordered_map<long int, bool> leftths = touchhs;
					inserthp(hpid, mink, node->left, leftths);
				}
				if (node->right->isPruned == false)
				{
					unordered_map<long int, bool> rightths = touchhs;
					inserthp(hpid, mink, node->right, rightths);
				}
				else if (node->left->isPruned == true && node->right->isPruned == true)
				{
					node->isPruned = true;
				}
			}
		}
	}
}


void cellTree::insert(cell *c, vector<long int> &hps, int mink, vector<cell*>& leaves)
{
	for (int pos = 0; pos < hps.size(); pos++)
	{	
	
		if (c->isPruned == false)
		{
			unordered_map<long int, bool> leftths;
			inserthp(hps[pos], mink, c, leftths);
		}
	}
	collectLeaf(c,leaves,mink);
}


void cellTree::collectLeaf(cell *c, vector<cell*>& leaves, const int& mink)
{
	leaves.clear();
	if (c->isPruned == false)
	{
		if (c->left != NULL )//&& root->left->isPruned ==  false)
		{
		//		cell* leaf = new cell();
			dfsTraversal(c->left, leaves);
		}
		if (c->right != NULL )//&& root->right->isPruned == false)
		{
		//	cell* leaf = new cell();
			dfsTraversal(c->right, leaves);
		}
	}
	delete c;
	//sort(leaves.begin(), leaves.end(), cellCompare());
}

void cellTree::dfsTraversal(cell* node, vector<cell*>& leaves)
{

	cell *new_leaf = new cell(node);

	if (node->left == NULL && node->right == NULL)
	{	
			leaves.push_back(new_leaf);	
	}
	else
	{
		if (node->left != NULL)
		{
			//cell* left = new cell();
			
			dfsTraversal(node->left, leaves);
			
		}
		if (node->right != NULL)
		{
			dfsTraversal(node->right, leaves);
		}
		//else
	//	{
	//		leaf->release();
	//		delete leaf;
	//	}
	}
	delete node;
	node = NULL;
}

void cellTree::opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult)
{
	for (int pos = 0; pos < hps.size(); pos++)
	{
		if (pos % 20 == 0 && pos != 0)
		{
			cout << pos << endl;
		}
		if (root->isPruned == false)
		{
			unordered_map<long int, bool> leftths;
			inserthp(hps[pos], mink, root->left, leftths);
			unordered_map<long int, bool> rightths;
			inserthp(hps[pos], mink, root->right, rightths);
		/*	cell* left = new cell();
			opt_inserthp(hps[pos], mink, root->right, left);
			cell* right = new cell();
			opt_inserthp(hps[pos], mink, root->right, right);*/
		}
	}
//	updateCellTree(root);
}

void cellTree::opt_inserthp(long int & hpid, const int mink, cell* node, cell* all)
{
	if (node->isPruned == false)
	{
		node->appendto(all);
		int indicator = isDominatorInserted(dagNode[hpid]->rDominator, all);
		if (indicator == -1) // found negative halfspace from hpid's dominators
		{
			node->AboveHP.insert(hpid);
			totalSpaceCost += 4;
			releaseCell(all);
			delete all;
		}
		else if (indicator == 1) // found positive halfspace from hpid's dominators
		{
			inserthp(hpid, mink, node, all->IntersectHP);
		}
		else if (indicator == 0) // go to its children
		{
			if (node->left->isPruned == false)
			{
				cell* left = new cell();
				all->appendto(left);
				opt_inserthp(hpid, mink, node->left, left);
			}
			if (node->right->isPruned == false)
			{
				opt_inserthp(hpid, mink, node->left, all);
			}
			else 
			{
				node->isPruned = true;
				all->release();
				delete all;
			}
		}
		else
		{
			cout << "there is a bug in opt_inserthp" << endl;
		}
	}
}

int cellTree::isDominatorInserted(unordered_set<long int>& rdominators, cell* all)
{
	for (auto iter = rdominators.begin(); iter != rdominators.end(); iter++)
	{
		if (all->AboveHP.find(*iter) != all->AboveHP.end())
		{
			return -1;
		}
		if (all->BelowHP.find(*iter) != all->BelowHP.end())
		{
			return 1;
		}
		if (all->IntersectHP.find(*iter) != all->IntersectHP.end())
		{
			return all->IntersectHP[*iter] == true ? 1 : -1;
		}
	}
	return 0;
}


void cellTree::updateCellTree(cell* node)
{
	if (node->isPruned == true || node == NULL)
		return ;
	else if (node->left == NULL&&node->right == NULL)
		return;
	else if (node->left->isPruned == true && node->right->isPruned == true)
	{
		node->isPruned = true;
		releaseCell(node->left);
		releaseCell(node->right);
		node->left = NULL;
		node->right = NULL;
		return;
	}
	else
	{
		updateCellTree(node->left);
		updateCellTree(node->right);
		if (node->left->isPruned == true && node->right->isPruned == true)
		{
			node->isPruned = true;
			releaseCell(node->left);
			releaseCell(node->right);
			node->left = NULL;
			node->right = NULL;
			return;
		}
	}
}




int cellTree::findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 1; 
	int isOptimal;
	cl.clear();
	cu.clear();

	// obtain min score for focal record
	for (int i = 0; i < dimen; i++)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = 1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
		//	cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cl.push_back(row[i]);
		}


		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = -1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
		//	cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cu.push_back(row[i]);
		}
	}
	return isOptimal;
}
