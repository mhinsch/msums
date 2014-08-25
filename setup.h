#ifndef SPMAIN_H
#define SPMAIN_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "sample.h"
#include "pairsample.h"
#include "anafunctors.h"
#include "groupana.h"

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


typedef Sample<string, char> StrSample;
typedef PairSample<string, char> PairStrSample;

class SSHandler : public AnalysisHandler<AnalysisBase<StrSample > >
	{
public:
	SSHandler();
	};

class PSHandler : public AnalysisHandler<AnalysisBase<PairStrSample > >
	{
public:
	PSHandler();
	};

template<class SAMPLEV>
class GroupStatHandler
	{
public:
	typedef vector<size_t> group_t;
	typedef AnalysisBase<SAMPLEV> analysis_t;

protected:
	vector<group_t> _groups;
	vector<string> _names;
	vector<analysis_t*> _stats;

public:
	void add(const string & name, const group_t & group)
		{
		analysis_t * ana = 0;

		if (name == "f3")
			ana = new PatF3<SAMPLEV>();
		else if (name == "f4")
			ana = new PatF4<SAMPLEV>();

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
