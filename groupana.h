#ifndef MULTIANA_H
#define MULTIANA_H

#include "anafunctors.h"
#include "stats_multi.h"

template<class SAMPLEV>
struct PatF3 : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples) const
		{
		VERIFY_MSG(samples.size() == 3, 
			"Error: Patterson's F3 requires 3 populations");

		return patterson_f3(
			samples[0]->alleles().begin(), samples[0]->alleles().end(),
			samples[1]->alleles().begin(), samples[1]->alleles().end(),
			samples[2]->alleles().begin(), samples[2]->alleles().end());
		}
	};


template<class SAMPLEV>
struct PatF4 : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples) const
		{
		VERIFY_MSG(samples.size() == 4, 
			"Error: Patterson's F4 requires 4 populations");

		return patterson_f4(
			samples[0]->alleles().begin(), samples[0]->alleles().end(),
			samples[1]->alleles().begin(), samples[1]->alleles().end(),
			samples[2]->alleles().begin(), samples[2]->alleles().end(),
			samples[3]->alleles().begin(), samples[3]->alleles().end());
		}
	};


template<class SAMPLEV>
struct FreqVarying : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples) const
		{
		if (samples.size() == 0)
			return 0;

		for (size_t i=0; i<samples.size(); i++)
			if (samples[i]->n_segregating_sites())
				return 1;

		for (size_t site=0; site<samples[0]->n_sites(); site++)
			for (size_t i=1; i<samples.size(); i++)
				if (samples[0]->alleles(site).find() != 
					samples[i]->alleles(site).find())
					return 1;
		
		return 0;		
		}
	};

#endif	// MULTIANA_H
