#ifndef SIDNA_H
#define SIDNA_H
 
#include <limits>
#include <vector>

using namespace std;


class Watterson
	{
protected:
	vector<float> _watt;

public:
	Watterson(int n)
		: _watt(n, 0.0)
		{}

	float operator()(int n)
		{
		float & w = _watt[n];

		if (!(w > 0))
			for (int i=1; i<n; i++)
				w += 1.0/i;

		return w;
		}
	};

inline double undefined()
	{
	return numeric_limits<double>::quiet_NaN();
	}

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

/*defines the structure used to store results from shared sites, fixed sites per locus*/
class SingleStats : public StatList<SStatsI, _sstats_last>
	{
protected:
	Watterson _watt;

public:
	SingleStats(int nLoci=0, int maxNSeq=0)
		: StatList<SStatsI, _sstats_last>(nLoci), _watt(maxNSeq) 		
		{ }

	void compute_pi_theta(const size_t locus,
		const vector<string> & seqlhs, int nspolyl);
	};


class PairStats : public StatList<PStatsI, _pstats_last>
	{
public:
	
	PairStats(int nLoci)
		: StatList<PStatsI, _pstats_last>(nLoci)
		{ }

	void compute_div_fst(const size_t locus,
		const vector<string> & seqAlhs, const vector<string> & seqBlhs,
		const SingleStats & statsA, const SingleStats & statsB,
		int nspolyl);

	void compute_polyl(const size_t locus,
		const vector<string> & seqAlhs, const vector<string> & seqBlhs, 
		int nspolyl, const string & seqOls, bool outgroup = false);

	};

float pearson_corr_pi(const SingleStats & statsA, const SingleStats & statsB);

#endif	// SIDNA_H
