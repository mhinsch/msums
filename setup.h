#ifndef SPMAIN_H
#define SPMAIN_H

#include <vector>
#include <map>
#include <set>
#include <string>

#include "anafunctors.h"
#include "sample.h"
#include "pairsample.h"

using namespace std;

template<class ANA>
struct AnalysisHandler
	{
public:
	typedef ANA analysis_t;

protected:
	map<string, analysis_t *> _stats;
	vector<string> _names;
	vector<analysis_t *> _active;
	vector<string> _active_names;
	
public:
	~AnalysisHandler()
		{
		// TODO: delete analysis objects
		}

	void add(const string & name, analysis_t * a)
		{
		_stats[name] = a;
		_names.push_back(name);
		}

	analysis_t * get(const string & name)
		{
		if (!exists(name))
			return 0;
		return _stats[name];
		}

	bool exists(const string & name) const
		{
		return _stats.count(name);
		}

	const vector<string> & names() const
		{ return _names; }

	void activate(set<string> & active)
		{
		_active.clear();
		_active_names.clear();

		for (size_t i=0; i<_names.size(); i++)
			if (active.count(_names[i]) && exists(_names[i]))
				{
				_active_names.push_back(_names[i]);
				_active.push_back(_stats[_names[i]]);
				}
		}

	const vector<analysis_t *> & active() const
		{ return _active; }

	const vector<string> & active_names() const
		{ return _active_names; }
	};


typedef Sample<string> StrSample;
typedef PairSample<string> PairStrSample;

class SSHandler : public AnalysisHandler<AnalysisBase<StrSample > >
	{
public:
	SSHandler()
		{
		typedef StrSample sample_t;

		typedef AnalysisWrapper<sample_t, double> SD;
		typedef AnalysisWrapper<sample_t, int> SI;

		add("pairdif", SI::create<&sample_t::sum_pairwise_differences>());
		add("segr", SI::create<&sample_t::n_segregating_sites>());
		add("singlet", SI::create<&sample_t::n_singletons>());
		add("thpi", SD::create<&sample_t::theta_pi>());
		add("thW", SD::create<&sample_t::theta_W>());
		add("flDstar", SD::create<&sample_t::fu_li_Dstar>());
		add("flFstar", SD::create<&sample_t::fu_li_Fstar>());
		add("tD", SD::create<&sample_t::tajima_D>());
		add("R2", SD::create<&sample_t::R2>());
		}
	};

class PSHandler : public AnalysisHandler<AnalysisBase<PairStrSample > >
	{
public:
	PSHandler()
		{
		typedef PairStrSample sample_t;

		typedef AnalysisWrapper<sample_t, double> SD;
		typedef AnalysisWrapper<sample_t, int> SI;

		add("d", SI::create<&sample_t::sum_pairwise_differences>());
		add("dn", SD::create<&sample_t::dn>());
		add("FST", SD::create<&sample_t::fst>());
		add("bialsites", SI::create<&sample_t::n_bial_sites>());
		add("multisites", SI::create<&sample_t::n_multi_sites>());
		add("sfA", SI::create<&sample_t::sfA>());
		add("sfB", SI::create<&sample_t::sfB>());
		add("sfout", SI::create<&sample_t::sfout>());
		add("sxA", SI::create<&sample_t::sxA>());
		add("sxB", SI::create<&sample_t::sxB>());
		add("sxAfB", SI::create<&sample_t::sxAfB>());
		add("sxBfA", SI::create<&sample_t::sxBfA>());
		add("ss", SI::create<&sample_t::ss>());
		add("Wald", SD::create<&sample_t::wald>());
		add("Rf", SI::create<&sample_t::Rf>());
		add("Rs", SI::create<&sample_t::Rs>());
		add("pattD", SD::create<&sample_t::patterson_D>());
		}
	};

template<class SAMPLEV>
class GroupStatHandler
	{
public:
	typedef vector<size_t> group_t;
	typedef AnalysisBase<SAMPLEI> analysis_t;

protected:
	vector<group_t> _groups;
	vector<string> _names;
	vector<analysis_t*> _stats;

public:
	void add(const string & name, const group_t & group)
		{
		analysis_t * ana = 0;

		if (name == "f3")
			ana = new PatF3<SAMPLEV>;
		else if (name == "f4")
			ana = new PatF4<SAMPLEV>;

		if (! ana)
			{
			// TODO: throw error
			return;
			}

		_stats.push_back(ana);
		_groups.push_back(group);
		_names.push_back(name);
		}

	const vector<string> & names() const
		{
		return _names;
		}

	const vector<group_t> & groups() const
		{
		return _groups;
		}

	const vector<analysis_t *> & stats() const
		{
		return _stats;
		}
	};


#endif	// SPMAIN_H
