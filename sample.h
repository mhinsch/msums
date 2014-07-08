#ifndef SAMPLE_H
#define SAMPLE_H

#include "stats.h"

template<class SEQ>
class Sample
	{
public:
	typedef SEQ sequence_t;
	typedef SEQ::value_type state_t;
	typedef std::pair<state_t, int> allele_count;
	typedef std::vector<allele_count> allele_set;

protected:
	// list of haplotypes in the sample
	vector<sequence_t> _sequences;
	// set of allele frequencies per site
	vector<allele_set> _alleles;
	// sites that are known to possibly be polymorphic
	size_t _n_poly_sites;


	template<class CONT>
	static CONT::iterator find_allele(state_t s, CONT & list)
		{
		CONT::iterator a=list.begin();

		for (; a!=list.end() && a->first != s; a++);

		return a;
		}
		
	int calc_sum_pairwise_differences() const
		{
		return sum_pairwise_differences(
			_alleles.begin(), _alleles.begin()+_n_poly_sites);
		}

	int calc_segregating_sites() const
		{
		return n_segregating_sites(_alleles.begin(), _alleles.begin()+_n_poly_sites);
		}

	int calc_n_singletons() const
		{
		return count_singletons(_alleles.begin(), _alleles.begin()+_n_poly_sites);
		}

	double calc_theta_pi() const
		{
		return theta_pi(size(), sum_pairwise_differences());
		}

	double theta_pi() const
		{
		return _theta_pi.ready() ? 
			_theta_pi :
			_theta_pi = theta_pi(size(), sum_pairwise_differences());
		}

	double calc_theta_W() const
		{
		return theta_W(size(), n_segregating_sites());
		}

	double calc_fu_li_Dstar() const
		{
		return fu_li_Dstar(size(), n_singletons(), n_segregating_sites());
		}

	double calc_fu_li_Fstar() const
		{
		return fu_li_Fstar(size(),n_singletons(), n_segregating_sites(), theta_pi());
		}

	double calc_tajima_D() const
		{
		return DTajima(size(), n_segregating_sites(), theta_pi());
		}


public:

	Sample()
		: sum_pairwise_differences(this), n_segregating_sites(this), 
		n_singletons(this), theta_pi(this), theta_W(this), 
		fu_li_Dstar(this), fu_li_Wstar(this), tajima_D(this)
		{}

	size_t n_sites() const
		{
		return _sequences[0].size();
		}
	size_t size() const
		{
		return _sequences.size();
		}

	void prepare_alleles_per_site()
		{
		_alleles.resize(_sequences[0].size());
		for (size_t i=0; i<_sequences.size(); i++)
			{
			const sequence_t & seq = _sequences[i];
			for (size_t site=0; site<_n_poly_sites; site++)
				{
				const state_t s = seq[i];
				allele_set::iterator a = find_allele(s, _alleles[i]);
				if (a != _alleles.end())
					a->second++;
				else
					_alleles.push_back(allele_count(s, 1));
				}
			}
		}

	bool n_alleles(size_t site)
		{
		return _alleles[site].size();
		}

	const allele_set & alleles(size_t site)
		{
		return _alleles[site];
		}

	const vector<allele_set> & alleles()
		{
		return _alleles;
		}

	typedef WrapCached<Sample> WC;

	mutable WC::CachedValue<int, &calc_sum_pairwise_differences>
		sum_pairwise_differences;

	mutable WC::CachedValue<int, &calc_segregating_sites>
		n_segregating_sites;

	mutable WC::CachedValue<int, &calc_n_singletons>
		n_singletons;

	mutable WC::CachedValue<double, &calc_theta_pi>
		theta_pi;
	
	mutable WC::CachedValue<double, &calc_theta_W>
		theta_W;

	mutable WC::CachedValue<double, &calc_fu_li_Dstar>
		fu_li_Dstar;

	mutable WC::CachedValue<double, &calc_fu_li_Wstar>
		fu_li_Wstar;

	mutable WC::CachedValue<double, &calc_tajima_D>
		tajima_D;
	};


#endif	// SAMPLE_H
