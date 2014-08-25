#ifndef SAMPLE_H
#define SAMPLE_H

#include <vector>

#include "stats_single.h"
#include "mathutil.h"
#include "util.h"
#include "alleleset.h"

using namespace std;


#define CACHED_OR_COMPUTE(name, calc) ( name.ready()?name:name(calc))

template<class SEQ, class STATE>
class Sample
	{
public:
	typedef SEQ sequence_t;
	typedef typename SEQ::value_type state_t;
	typedef AlleleSet<STATE> allele_set;

protected:
	// list of haplotypes in the sample
	vector<sequence_t> _sequences;
	// set of allele frequencies per site
	vector<allele_set> _alleles;

	size_t _tot_n_sites;

	bool _prepared;

	typedef CachedValue<double> CVD;
	typedef CachedValue<int> CVI;

	mutable CVI _spwd, _nss, _nst;
	mutable CVD _theta_pi, _theta_W, _flD, _flF, _tajD, _R2;
	mutable CachedValue<vector<int> > _nstps;

	template<class CONT>
	static typename CONT::iterator find_allele(state_t s, CONT & list)
		{
		typename CONT::iterator a=list.begin();

		for (; a!=list.end() && a->first != s; a++);

		return a;
		}
		
public:

	Sample()
		: _sequences(0), _alleles(0), _prepared(false)
		{
		}


	Sample(size_t n_seq, const sequence_t & ini)
		: _sequences(n_seq, ini), _alleles(0), _prepared(false)
		{
		}

	void reset()
		{
		//cerr << "reset\n";
		// reset all analysis objects
		_spwd.reset();
		_nss.reset();
		_nst.reset();
		_nstps.reset();
		_theta_pi.reset();
		_theta_W.reset();
		_flD.reset();
		_flF.reset();
		_tajD.reset();
		_R2.reset();

		// clear sequence data
		for (size_t sq=0; sq<_sequences.size(); sq++)
			_sequences[sq].clear();
		// clear allele counts
		for (size_t st=0; st<_alleles.size(); st++)
			_alleles[st].clear();

		_prepared = false;
		}

	size_t size() const
		{ return _sequences.size(); }
	void resize(size_t s, const sequence_t & ini)
		{
		_sequences.resize(s, ini);
		}

	sequence_t & sequence(size_t i)
		{ return _sequences[i]; }
	const sequence_t & sequence(size_t i) const
		{ return _sequences[i]; }

	size_t n_sites() const
		{ return _sequences[0].size(); }


	size_t tot_n_sites() const
		{
		return _tot_n_sites;
		}
	void set_tot_n_sites(size_t n)
		{
		_tot_n_sites = n;
		}

	void prepare_alleles_per_site()
		{
		// let's not do this twice
		if (_prepared) return;
		_prepared = true;

		// one alleleset per site
		_alleles.resize(_sequences[0].size(), allele_set());

		// record alleles for each sequence
		for (size_t sq=0; sq<_sequences.size(); sq++)
			for (size_t site=0; site<_sequences[sq].size(); site++)
				_alleles[site].add(_sequences[sq][site]);

		// now count alleles
		for (size_t s=0; s<_alleles.size(); s++)
			_alleles[s].do_count();
		}

	bool n_alleles(size_t site) const
		{ return _alleles[site].count(); }

	const allele_set & alleles(size_t site) const
		{ return _alleles[site]; }

	const vector<allele_set> & alleles() const
		{
		return _alleles;
		}

	int sum_pairwise_differences() const
		{
		return CACHED_OR_COMPUTE(_spwd, pairwise_differences(
				_alleles.begin(), _alleles.end()));
		}

	int n_segregating_sites() const
		{
		return CACHED_OR_COMPUTE(_nss, ::n_segregating_sites(
				_alleles.begin(), _alleles.end()));
		}

	int n_singletons() const
		{
		return CACHED_OR_COMPUTE(_nst, count_singletons(
				_alleles.begin(), _alleles.end()));
		}

	const vector<int> & n_singletons_by_seq() const
		{
		if (_nstps.ready())
			return _nstps;

		vector<int> & n = _nstps.get();
		n.resize(_sequences.size());
		for (int i=0; i<n.size(); i++)
			n[i] = count_singletons_for_seq(
				_sequences[i].begin(), _sequences[i].end(), _alleles.begin());

		_nstps.set_ready();
		
		return _nstps;
		}

	double theta_pi_s() const
		{
		return theta_pi() / _tot_n_sites;
		}

	double theta_pi() const
		{
		return CACHED_OR_COMPUTE(_theta_pi, ::theta_pi(
				size(), sum_pairwise_differences()));
		}

	double theta_W_s() const
		{
		return theta_W() / _tot_n_sites;
		}

	double theta_W() const
		{
		return CACHED_OR_COMPUTE(_theta_W, ::theta_W(
				size(), n_segregating_sites()));
		}

	double fu_li_Dstar() const
		{
		return CACHED_OR_COMPUTE(_flD, ::fu_li_Dstar(
				size(), n_singletons(), n_segregating_sites()));
		}

	double fu_li_Fstar() const
		{
		return CACHED_OR_COMPUTE(_flF, ::fu_li_Fstar(
				size(), n_singletons(), n_segregating_sites(), theta_pi()));
		}

	double tajima_D() const
		{
		return CACHED_OR_COMPUTE(_tajD, DTajima(
				size(), n_segregating_sites(), theta_pi()));
		}

	double R2() const
		{
		const vector<int> & ns = n_singletons_by_seq();

		return CACHED_OR_COMPUTE(_R2, ::R2(
				ns.begin(), ns.end(), theta_pi(), n_segregating_sites()));
		}
	};


#endif	// SAMPLE_H
