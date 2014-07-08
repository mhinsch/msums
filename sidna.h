#ifndef SIDNA_H
#define SIDNA_H
 
#include <vector>
#include <string>


#include "util.h"

using namespace std;


template<typename T>
class Stats
	{
protected:
	vector<T> _values;
	T _std;
	T _mean;

public:
	Stats(size_t num=0)
		{
		_values.resize(num, undefined());
		}
	
	void resize(size_t s)
		{
		_values.resize(s, undefined());
		}

	void reset()
		{
		fill(_values.begin(), _values.end(), undefined());
		}

	void compute_aggr()
		{
		T sqr_sum = 0;
		size_t n = 0;
		_mean = 0;

		for (size_t i=0; i<_values.size(); i++)
			{
			// NAN
			if (_values[i] != _values[i])
				continue;
			n++;
			_mean += _values[i];
			sqr_sum += _values[i] * _values[i];
			}

		if (n<2)
			return;
		
		const T ssq = sqr_sum - _mean*_mean/n;
		_std = sqrt(ssq / T(n-1));

		_mean /= n;
		}

	T std() const
		{
		return _std;
		}
	T mean() const
		{
		return _mean;
		}

	T & operator[](size_t i)
		{
		return _values[i];
		}

	const T & operator[](size_t i) const
		{
		return _values[i];
		}

	size_t size() const
		{
		return _values.size();
		}
	};

typedef Stats<double> RStats;

enum SStatsI {
	_sumpairdif,
	_theta,
	_D,
	_pi,
	_Dstar,
	_Fstar,
	_sstats_last};

enum PStatsI {
	_d,
	_dn,
	_FST,
	_bialsites,
	_multisites,
	_sfA,
	_sfB,
	_sfout,
	_sxA,
	_sxB,
	_sxAfB,
	_sxBfA,
	_ss,
	_Wald,
	_pstats_last};


inline bool requiresSinglePopResults(PStatsI idx)
	{
	return idx <= _FST;
	}


template<typename IDX, IDX MAX>
class StatList
	{
protected:
	RStats _stats[MAX];

public:
	typedef IDX idx_type;
	
	StatList(int nLoci)
		{
		for (int i=0; i<MAX; i++)
			_stats[i] = RStats(nLoci);
		}

	void reset()
		{
		for (int i=0; i<MAX; i++)
			_stats[i].reset();
		}

	void set(IDX idx, int locus, double value)
		{
		_stats[idx][locus] = value;
		}

	double get(IDX idx, int locus) const
		{
		return _stats[idx][locus];
		}

	void compute_aggr()
		{
		for (int idx=0; idx<int(MAX); idx++)
			_stats[IDX(idx)].compute_aggr();
		}

	RStats & operator[](IDX idx)
		{
		return _stats[idx];
		}

	const RStats & operator[](IDX idx) const
		{
		return _stats[idx];
		}
	};


double pearson_corr_pi(const SingleStats & statsA, const SingleStats & statsB);

#endif	// SIDNA_H
