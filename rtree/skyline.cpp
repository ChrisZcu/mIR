#include "skyline.h"
typedef multimap<float, long int>::value_type PfltINT;
typedef multimap<long int, VirtualRNode*>::value_type PintVRN;

extern unordered_map<long int, RtreeNode*> ramTree;

void rtreeRAM(Rtree& rt, unordered_map<long, RtreeNode*>& ramTree)
{
	ramTree.clear();
	queue<long> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long pageID;
	while (!H.empty())
	{
		pageID = H.front();
		
		H.pop();
		node = rt.m_memory.loadPage(pageID);
		ramTree[pageID] = node;
		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
	}
}

int countkDominator(const int dimen, const float pt[], vector<long> kskyband, float* PG[])
{
	vector<long>::iterator iter;
	if (kskyband.size() == 0)
		return false;

	int count = 0;
	for (iter = kskyband.begin(); iter != kskyband.end(); iter++)
	{
		long pid = *iter;
		float s[MAXDIMEN];
		bool dominated = true;
		for (int i = 0; i < dimen; i++)
		{
			if (PG[pid][i] + SIDELEN < pt[i])
			{
				dominated = false;
				break;
			}
		}
		if (dominated)
			count++;
	}
	return count;
}


float minDist(float p1[], float p2[], int dimen)
{
	float mindist = 0;

	for (int i = 0; i < dimen; i++)
	{
		float dist = p1[i] - p2[i];
		mindist += (dist * dist);
	}
	return (float)sqrt(mindist);
}




void GetSkylines(const int dimen, Rtree& a_rtree, std::multimap<long int, VirtualRNode*>& NonResultEntry, std::vector<long int>& PrunedNodes, set<long int>& a_skylines, float* PG[])
{
	multimap<float, long int> H0;
	vector<long int> skylines;

	float pt[MAXDIMEN];
	bool isAnObject;
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	RtreeNodeEntry* e0;

	set<long int>::iterator s0Iter, s1Iter;
	multimap<long int, VirtualRNode*>::iterator nretIter;
	vector<long int>::iterator pIter;

	for (pIter = PrunedNodes.begin(); pIter != PrunedNodes.end(); pIter++)
	{
		RtreeNode* rNode = a_rtree.m_memory.loadPage(*pIter);
		e0 = rNode->genNodeEntry();
		for (int i = 0; i < dimen; i++)
		{
			pt[i] = e0->m_hc.getUpper()[i];
		}
		mindist = minDist(pt, ORIGNIN, dimen);
		H0.insert(PfltINT(mindist, *pIter));
		delete rNode;
		delete e0;
	}
	PrunedNodes.clear();

	vector<long int> nIDToDelete;
	for (nretIter = NonResultEntry.begin(); nretIter != NonResultEntry.end(); nretIter++)
	{
		if (a_skylines.size() > 0)
		{
			s0Iter = a_skylines.find(nretIter->first - MAXPAGEID);
			if (s0Iter != a_skylines.end())
			{
				nIDToDelete.push_back(nretIter->first);
				continue;
			}
		}

		for (int i = 0; i < dimen; i++)
		{
			pt[i] = (nretIter->second)->m_entry[0]->m_hc.getUpper()[i];
		}
		mindist = minDist(pt, ORIGNIN, dimen);
		H0.insert(PfltINT(mindist, nretIter->first));
	}

	for (pIter = nIDToDelete.begin(); pIter != nIDToDelete.end(); pIter++)
	{
		nretIter = NonResultEntry.find(*pIter);
		delete nretIter->second;
		NonResultEntry.erase(nretIter);
	}

	multimap<float, long int>::iterator fIter;
	while (H0.size() != 0)
	{
		isAnObject = false;

		fIter = H0.begin();
		long int pageid = fIter->second;
		float distTemp = fIter->first;
		H0.erase(fIter);

		VirtualRNode* VirNode = new VirtualRNode;
		RtreeNodeEntry* e0;            // create node entry e0 for node n, so that its MBR can be obtained    	
		if (pageid >= MAXPAGEID)     //current element in H is a data entry node (for distinction,its pageid is equal to m_id+MAXPAGED)
		{
			isAnObject = true;
			nretIter = NonResultEntry.find(pageid);
			if (nretIter == NonResultEntry.end())
			{
				cout << "Error! there is no node " << pageid << " in NonResultEntry!" << endl;
			}
			else
			{
				VirNode->copyData(nretIter->second);
				if (VirNode->m_usedspace > 1)
				{
					cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
				}
				else
					e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node           	
			}
			pageid = pageid - MAXPAGEID;
		}
		else
		{
			RtreeNode* node = a_rtree.m_memory.loadPage(pageid);
			VirNode->copyData(*node);
			VirNode->copyEntries(*node, node->m_usedspace);
			e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node        
			delete node;
		}

		if (isAnObject)
		{
			for (int i = 0; i < dimen; i++)
				pt[i] = VirNode->m_entry[0]->m_hc.getLower()[i] + SIDELEN;
		}
		else
		{
			for (int i = 0; i < dimen; i++)
				pt[i] = e0->m_hc.getUpper()[i];
		}
		bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

		if (!dominated)
		{
			if (VirNode->isLeaf())
			{
				if (VirNode->m_usedspace>1)
				{
					for (int i = 0; i < VirNode->m_usedspace; i++)
					{
						for (int j = 0; j < dimen; j++)
						{
							pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];
						}

						bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

						long int NewPageId = VirNode->m_entry[i]->m_id + MAXPAGEID;
						VirtualRNode* node = new VirtualRNode;
						RtreeNodeEntry* newEntry = VirNode->m_entry[i]->clone();
						node->insertEntry(newEntry);
						NonResultEntry.insert(PintVRN(NewPageId, node));
						delete newEntry;

						if (dominated)
							continue;

						mindist = minDist(pt, ORIGNIN, dimen);
						H0.insert(PfltINT(mindist, NewPageId));
					}
				}
				else
				{
					skylines.push_back(pageid);
				}
			}
			else
			{
				for (int i = 0; i < VirNode->m_usedspace; i++)
				{
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];
					}

					bool dominated = IsDominatedBy(dimen, pt, skylines, PG);
					if (!dominated)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						H0.insert(PfltINT(mindist, VirNode->m_entry[0]->m_id));
					}
					else
					{
						PrunedNodes.push_back(VirNode->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			if (VirNode->isLeaf())
			{
				if (VirNode->m_usedspace > 1)
				{
					PrunedNodes.push_back(pageid);
				}
			}
			else
			{
				PrunedNodes.push_back(pageid);
			}
		}
		delete VirNode;
		delete e0;
	}
	a_skylines.clear();
	for (vector<long>::iterator iter = skylines.begin(); iter != skylines.end(); iter++)
	{
		a_skylines.insert(*iter);
	}
}

bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[])
{
	vector<long>::iterator iter;
	if (a_skylines.size() == 0)
		return false;

	for (iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		long pid = *iter;
		float s[MAXDIMEN];
		bool dominated = true;
		for (int i = 0; i < dimen; i++)
		{
			if (PG[pid][i] + SIDELEN < pt[i])
			{
				dominated = false;
				break;
			}
		}
		if (dominated)
			return dominated;
	}
	return false;
}

bool isDominateByFocal(const int dimen, const float pt[], Point& focal)
{
	bool dominated = true;
	for (int i = 0; i < dimen; i++)
	{
		if (focal.m_coor[i] + SIDELEN < pt[i])
		{
			dominated = false;
			break;
		}
	}
	return dominated;
}

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<long int, RtreeNodeEntry*> RecordEntry;
	multimap<long int, RtreeNodeEntry*>::iterator recordIter;
	multimap<float, int> heap;
	int NegPageid;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	node = a_rtree.m_memory.loadPage(a_rtree.m_memory.m_rootPageID);

	if (node->isLeaf())
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
			
			NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
			RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
			RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
		}
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
		}
	}
	delete node;

	while (heap.size() != 0)
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			recordIter = RecordEntry.find(pageID);
			if (recordIter != RecordEntry.end())
			{
				for (int d = 0; d < dimen; d++)
				{
					pt[d] = recordIter->second->m_hc.getLower()[d] + SIDELEN;
				}
				if (countkDominator(dimen, pt, kskyband, PG) <= k && isDominateByFocal(dimen, pt, a_pt) == false)
				{
					kskyband.push_back(pageID - MAXPAGEID);
				}
			}
			else
			{
				cout << "an error incured in Getkskyband" << endl;
				exit(0);
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
					}
					if (countkDominator(dimen, pt, kskyband, PG) <= k) 
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
						
						NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
						RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
						RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
			delete node;
		}
	}
}

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<float, int> heap;
	multimap<float, int>::iterator heapIter;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;

	heap.insert(PfltINT(INFINITY, a_rtree.m_memory.m_rootPageID));

	while (!heap.empty())
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			for (int d = 0; d < dimen; d++)
				pt[d] = (PG[pageID - MAXPAGEID][d] + PG[pageID - MAXPAGEID][d + dimen])/2;
			if (countkDominator(dimen, pt, kskyband, PG) <= k)
			{
				kskyband.push_back(pageID - MAXPAGEID);
			}
		}
		else
		{
			node = ramTree[pageID];
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
					}
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
		}
	}
}

void printSkyBand(float* PointSet[], vector<long int>& skyband, int& dim, char* filename) {
   ofstream fout(filename);
  if (!fout.is_open()) {
    std::cerr << "in data_generator.printVVF()";
    std::cerr << filename << " not open" << std::endl;
    exit(-1);
  }
  
  fout.precision(4);
  for(int i=0; i<skyband.size(); i++) {
    fout << i+1 << " ";
    for(int j=0; j<2*dim; j++) {
      fout << PointSet[skyband[i]][j] << " ";
    }
    fout << endl;
  }
}