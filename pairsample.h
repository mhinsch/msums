#ifndef GROUPSAMPLE_H
#define GROUPSAMPLE_H

#include <vector>

#include "sample.h"
#include "stats_multi.h"

template<class SEQ>
class PairSample
	{
public:
	typedef Sample<SEQ> sample_t;

protected:
	const sample_t & _p1, & _p2;

	mutable CachedValue<int> _spd;
	mutable CachedValue<double> _fst;
	mutable CachedValue<double> _patterson_D;
	mutable CachedValue<PolyL> _poly_l;
	mutable CachedValue<pair<int, int> > _navascues_R;
	mutable CachedValue<pair<double, double> > _navascues_W;

public:
	PairSample(const sample_t & p1, const sample_t & p2)
		: _p1(p1), _p2(p2)
		{
		}

	int sum_pairwise_differences() const
		{
		if (_spd.ready())
			return _spd;

		int n_poly_sites = std::max(_p1.n_sites(), _p2.n_sites());
		return _spd = sum_pair_diff(
			_p1.alleles().begin(), _p1.alleles().begin() + n_poly_sites,
			_p2.alleles().begin(), _p2.alleles().begin() + n_poly_sites);
		}

	double fst() const
		{
		if (_fst.ready())
			return _fst;

		const int n_pairs = _p1.size() * _p2.size();
		const int count_pair_piT = 
			_p1.size() * (_p1.size() - 1) / 2 + _p2.size() * (_p2.size() - 1) / 2 + 
			n_pairs;

		const double piT = 
			(_p1.sum_pairwise_differences() + _p2.sum_pairwise_differences() + 
			sum_pairwise_differences()) / (double) count_pair_piT;

		return _fst(piT >= 1.0e-7 ?
			1.0 - (_p1.theta_pi()+_p2.theta_pi())/2.0 / piT : 
			undefined());
		}

	double patterson_D() const
		{
		return CACHED_OR_COMPUTE(_patterson_D, ::patterson_D(
				_p1.alleles().begin(), _p1.alleles().end(),
				_p2.alleles().begin(), _p2.alleles().end()));
		}

	// name?
	double dn() const
		{
		return sum_pairwise_differences() - (_p1.theta_pi() + _p2.theta_pi())/2;
		}

	const PolyL & poly_l() const
		{
		if (!_poly_l.ready())
			{
			_poly_l.get().compute(
				_p1.alleles().begin(), _p1.alleles().end(),
				_p2.alleles().begin(), _p2.alleles().end());
			_poly_l.set_ready();
			}

		return _poly_l;
		}

	int n_bial_sites() const
		{ return poly_l().bialsites; }

	int n_multi_sites() const
		{ return poly_l().multisites; }

	int sfout() const
		{ return poly_l().sfout; }

	int sfA() const
		{ return poly_l().sfA; }

	int sfB() const
		{ return poly_l().sfB; }

	int sxA() const
		{ return poly_l().sxA; }

	int sxB() const
		{ return poly_l().sxB; }

	int sxAfB() const
		{ return poly_l().sxAfB; }

	int sxBfA() const
		{ return poly_l().sxBfA; }

	int ss() const
		{ return poly_l().ss; }

	const pair<int, int> & navascues_R() const
		{
		return CACHED_OR_COMPUTE(_navascues_R, ::navascues_R(
				_p1.alleles().begin(), _p1.alleles().end(),
				_p2.alleles().begin(), _p2.alleles().end()));
		}

	int Rf() const
		{ return navascues_R().first; }

	int Rs() const 
		{ return navascues_R().second; }

	const pair<double, double> & navascues_W() const
		{
		return CACHED_OR_COMPUTE(_navascues_W, navascues_W_01(
				_p1.alleles().begin(), _p1.alleles().end(),
				_p2.alleles().begin(), _p2.alleles().end()));
		}

	double Wx2s1() const
		{ return navascues_W().first; }

	double Wx1s2() const 
		{ return navascues_W().second; }
	};


#endif	// GROUPSAMPLE_H
