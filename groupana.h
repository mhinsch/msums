#ifndef MULTIANA_H
#define MULTIANA_H

#include "anafunctors.h"
#include "stats_multi.h"

template<class SAMPLEV>
struct PatF3 : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples)
		{
		return patterson_f3(
			samples[0].alleles().begin(), samples[0].alleles().end(),
			samples[1].alleles().begin(), samples[1].alleles().end(),
			samples[2].alleles().begin(), samples[2].alleles().end());
		}
	};


template<class SAMPLEV>
struct PatF4 : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples)
		{
		return patterson_f4(
			samples[0].alleles().begin(), samples[0].alleles().end(),
			samples[1].alleles().begin(), samples[1].alleles().end(),
			samples[2].alleles().begin(), samples[2].alleles().end(),
			samples[3].alleles().begin(), samples[3].alleles().end());
		}
	};

#endif	// MULTIANA_H
