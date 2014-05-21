#ifndef SPMAIN_H
#define SPMAIN_H


#include <cmath>

#include <vector>
#include <map>

#include "sidna.h"

using namespace std;


struct StatDict
	{
	vector<string> singles;
	vector<string> pairs;
	map<string, SStatsI> idxSingles;
	map<string, PStatsI> idxPairs;

	StatDict()
		{
		addS("sumpairdif");
		addS("theta");
		addS("D");
		addS("pi");

		addP("d");
		addP("dn");
		addP("FST");
		addP("bialsites");
		addP("multisites");
		addP("sfA");
		addP("sfB");
		addP("sfout");
		addP("sxA");
		addP("sxB");
		addP("sxAfB");
		addP("sxBfA");
		addP("ss");
		addP("Wald");
		}

	void addS(const string & name)
		{
		idxSingles[name] = SStatsI(singles.size());
		singles.push_back(name);
		}

	void addP(const string & name)
		{
		idxPairs[name] = PStatsI(pairs.size());
		pairs.push_back(name);
		}
	};


#endif	// SPMAIN_H
