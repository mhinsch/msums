#ifndef SPMAIN_H
#define SPMAIN_H


#include <cmath>

#include <vector>
#include <map>

#include "sample.h"
#include "anafunctors.h"

using namespace std;


struct AnaHandler
	{
	map<string, AnalysisBase *> single_stats;
	vector<string> ss_names;
	vector<AnalysisBase *> use_ss;
	vector<

	StatDict()
		{
		add_ss("pairdif", create_analysis(&Sample::sum_pairwise_differences));
		add_ss("segr", create_analysis(&Sample::n_segregating_sites));
		add_ss("singlet", create_analysis(&Sample::n_singletons));
		add_ss("thpi", create_analysis(&Sample::theta_pi));
		add_ss("thW", create_analysis(&Sample::theta_W));
		add_ss("flDstar", create_analysis(&Sample::fu_li_Dstar));
		add_ss("flWstar", create_analysis(&Sample::fu_li_Wstar));
		add_ss("tD", create_analysis(&Sample::tajima_D));


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
	
	void add_ss(const string & name, AnalysisBase * a)
		{
		single_stats[name] = a;
		ss_names.push_back(name);
		}

	bool exists_ss(const string & name) const
		{
		return single_stats.has_key(name);
		}

	void create_ss_lists(const set<string> & names)
		{


	};


#endif	// SPMAIN_H
