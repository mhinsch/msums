#include "setup.h"


SSHandler::SSHandler()
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

PSHandler::PSHandler()
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
	add("Wx2s1", SD::create<&sample_t::Wx2s1>());
	add("Wx1s2", SD::create<&sample_t::Wx1s2>());
	add("pattD", SD::create<&sample_t::patterson_D>());
	}



