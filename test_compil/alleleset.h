#ifndef ALLELESET_H
#define ALLELESET_H

#include <vector>

// Default implementation for AlleleSet. Note that this will become rather
// inefficient for larger alphabets (e.g. microsatellites).
template<class STATE>
class AlleleSet
	{
public:
	typedef STATE state_t;

protected:
	int _n_alleles;
	std::vector<int> _count;

public:
	AlleleSet()
		: _n_alleles(0), _count(0)
		{}

	size_t s2i(const state_t & a) const
		{
		return a;
		}

	state_t i2s(size_t i) const
		{
		return i;
		}

	void add(const state_t & a)
		{
//		cerr << "ADD " << a << ";";
		const size_t i = s2i(a);
		if (_count.size() <= i)
			_count.resize(i+1, 0);
		_count[i]++;
		}

	int count(const state_t & a) const
		{
		const size_t i = s2i(a);
		return _count.size() > i ? _count[i] : 0;
		}

	// unsafe version
	int operator[](const state_t & a) const
		{
		return _count[s2i(a)];
		}

	// find first allele with count>0
	state_t find(size_t start = 0) const
		{
		for (size_t i=start; i<_count.size(); i++)
			if (_count[i] != 0)
				return i2s(i);
		
		return i2s(_count.size());
		}

	size_t size() const
		{
		return _count.size();
		}

	void resize(size_t size)
		{
		_count.resize(size);
		}

	void clear()
		{
		_count.clear();
		_n_alleles = 0;
		}

	void do_count()
		{
		_n_alleles = 0;

		for (size_t i=0; i<_count.size(); i++)
			if (_count[i] != 0)
				_n_alleles++;

		//cout << "na: " << _n_alleles << endl;
		}

	int product(const AlleleSet<STATE> & o) const
		{
		int p = 0;

		for (size_t i=0; i<_count.size() && i<o._count.size(); i++)
			p += _count[i] * o._count[i];
		
		return p;
		}

	int pairwise_difference() const
		{
		int diff = 0;

		for (size_t i=0; (i+1)<_count.size(); i++)
			for (size_t j=i+1; j<_count.size(); j++)
				diff += _count[i] * _count[j];

		//cerr << "diff: " << diff << endl;		

		return diff;
		}

	int n_singletons() const
		{
		int n = 0;

		for (size_t i=0; i<_count.size(); i++)
			if (_count[i] == 1) n++;

		return n;
		}

	bool is_singleton(const state_t & s) const
		{
		return count(s) == 1;
		}

	int n_alleles() const
		{
		return _n_alleles;
		}
	};

#endif	// ALLELESET_H
