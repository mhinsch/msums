#ifndef GROUPSAMPLE_H
#define GROUPSAMPLE_H

#include <vector>

#include "sample.h"
#include "stats_multi.h"


class PairSample
	{
protected:
	const Sample & _p1, & _p2;

	int calc_sum_pair_diff() const
		{
		int n_poly_sites = std::max(_p1.n_poly_sites(), _p2.n_poly_sites());
		return sum_pair_diff(
			_p1.alleles().begin(), _p1.alleles().begin() + n_poly_sites,
			_p2.alleles().begin(), _p2.alleles().begin() + n_poly_sites)
		}

	double calc_fst() const
		{
		const int n_pairs = _p1.size() * _p2.size();
		const int count_pair_piT = 
			_p1.size() * (_p1.size() - 1) / 2 + _p2.size() * (_p2.size() - 1) / 2 + 
			n_pairs;

		const double piT = 
			(_p1.sum_pair_diff() + _p2.sum_pair_diff() + 
			sum_pair_diff()) / (double) count_pair_piT;

		return piT >= 1.0e-7 ?
			1.0 - (_p1.theta_pi()+_p2.theta_pi())/2.0 / piT : 
			undefined();
		}
public:
	PairSample(const Sample & p1, const Sample & p2)
		: _p1(p1), _p2(p2)
		{
		}

	// name?
	double dn() const
		{
		return sum_pairwise_differences() - (_p1.theta_pi() + _p2.theta_pi())/2;
		}

	typedef WrapCached<PairSample> WC;

	mutable WC::CachedValue<int, &calc_sum_pair_diff>
		sum_pairwise_differences;

	mutable WC::CachedValue<double, &calc_fst>
		fst;
	};



#endif	// GROUPSAMPLE_H
