#include <cstdlib>

#include <vector>
#include <iostream>

#include "spmain.h"
#include "sidna.h"

using namespace std;


/** Checks whether number of alleles at @site is > 2. */
bool check_multi(const Sample & sequences, 
	size_t site, biallelic_sites & result)
	{
	const char first_base = sequences[0][site];

	result.n = 1;
	result.b1 = first_base;

	for (size_t h=1; h<sequences.size(); h++) 
		{	/*starts loop over haplotypes*/
		char current_site = sequences[h][site];

		if (current_site != first_base) 
			{	/* if base different than first haplotype*/
			if (result.n == 1) 
				{		/* the site was not known as biallelic within A*/
				result.n = 2;
				result.b2 = current_site;
				} 
			else if (result.n == 2) 
				{		/* the site was known as biallelic*/
				if (current_site != result.b2) 
					{	/*more than 2 different bases then multiallelic*/
					result.n = 3;
					return true;
					}
				}	/*end else if*/
			}	/*end if current_site != first_base*/
		}	/*end loop over haplotypes at species A*/
	
	return false;
	}




typedef vector<char> wald_type;
double WaldWolfowitz(const wald_type & IOseq, int stats);




void PairStats::compute_div_fst(const size_t locus,
	const Sample & seqAlhs, const Sample & seqBlhs,
	const SingleStats & statsA, const SingleStats & statsB,
	int nspolyl)
	{
	const int nseqA = seqAlhs.size();
	const int nseqB = seqBlhs.size();
	const int nSites = seqAlhs[0].size();

	/*now computes average pairwise divergence between species A and B*/
	const int count_pairAB = nseqA * nseqB;
	const double pair_div = avg_pair_div(seqAlhs, seqBlhs, nspolyl) / double(nSites);

	set(_d, locus, pair_div / double(count_pairAB));
	set(_dn, locus, 
		get(_d, locus) - (statsA.get(_pi, locus) + statsB.get(_pi, locus))/2.0);

	const int count_pair_piT = 
		nseqA * (nseqA - 1) / 2 + nseqB * (nseqB - 1) / 2 + count_pairAB;
	
	const double piT = 
		(statsA.get(_sumpairdif, locus) + statsB.get(_sumpairdif, locus) + 
		pair_div) / (double) count_pair_piT;

	if (piT >= 1.0e-7) 
		set(_FST, locus, 
			(piT - (statsA.get(_pi, locus)+statsB.get(_pi, locus))/2.0)/piT);
	}

void PairStats::compute_polyl(const size_t locus,
	const Sample & seqAlhs, const Sample & seqBlhs, 
	int nspolyl, const Sequence & seqOls, bool outgroup)
	{
	const int nSites = seqAlhs[0].size();

	int count_bialsites=0;
	int count_multisites=0;
	int count_sfA=0;
	int count_sfB=0;
	int count_sfout=0;
	int count_sxA=0;
	int count_sxB=0;
	int count_sxAfB=0;
	int count_sxBfA=0;
	int count_ss=0;

	/*vector of runs of fixed and polymorphic sites according to Miguel Navascues*/
	wald_type Walds(nSites);

	for (int s=0; s<nspolyl; s++) 
		{	/*starts loop over polymorphic sites at locus l*/
		biallelic_sites biA, biB;

		const bool is_multiallelic = 
			check_multi(seqAlhs, s, biA) || check_multi(seqBlhs, s, biB);

		if (is_multiallelic) 
			{	/*multiple hits*/
			count_multisites++;
			continue;	/*skip to next site*/
			}

		// only fixed or biallelic arrives here

		const char outgroup = seqOls[s];

		// now test whether multiple hits and if not 
		// determine the type of polymorphism
		const int code = biA.n + biB.n*2;

		switch (code)
			{
		case 3:	// a==1, b==1
			// all identical
			if (biA.b1==biB.b1 && biA.b1==outgroup)
				continue;
			// all different => multiallelic
			if (biA.b1!=biB.b1 && biA.b1!=outgroup && biB.b1!=outgroup)
				break;

			// two alleles
			count_bialsites++;

			if (biA.b1 == biB.b1)
				/*site fixed along outgroup*/
				count_sfout++;
			else	// exactly one of a,b identical to outgroup
				{	
				Walds[s] = true;
				biA.b1 == outgroup ?
					count_sfB++ :	// base in B is derived
					count_sfA++; 	// base in A is derived
				}

			continue;
		case 4: // a==2, b==1
			if ( (biB.b1==biA.b1 || biB.b1==biA.b2) &&
				(outgroup==biB.b1 || outgroup==biA.b1 || outgroup==biA.b2) ) 
				{	/* fixed in B for one of the two bases in A*/
				Walds[s] = false;
				count_bialsites++;

				if (outgroup == biB.b1) 	/*unique derived polymorphism in A*/
					count_sxA++;
				else					/*unique ancestral polymorphism in A*/
					count_sxAfB++;

				continue;
				} 
			break;
		case 5: // a==1, b==2
			// fixed in A for one out of B and one of them equal to outgroup
			if ( (biA.b1==biB.b1 || biA.b1==biB.b2) &&
				(outgroup==biA.b1 || outgroup==biB.b1 || outgroup==biB.b2) )
				{
				Walds[s] = false;
				count_bialsites++;

				if (outgroup == biA.b1)	// unique derived polymorphism in B
					count_sxB++;
				else 					// unique ancestral polymorphism in B
					count_sxBfA++;

				continue;
				} 
			// base fixed in A is different from both in B
			// or outgroup different from all
			break;
		case 6:	// a==2, b==2
			// shared polymorphism and outgroup present
			if (  ( (biA.b1==biB.b1 && biA.b2==biB.b2) ||
					(biA.b1==biB.b2 && biA.b2==biB.b1) )  &&
				(outgroup==biA.b1 || outgroup==biA.b2)  )
				{
				count_ss++;
				Walds[s] = false;
				count_bialsites++;
				continue;
				} 
			/* polymorphism for different pairs of bases in A and B*/
			break;
		default:
			cerr << "Error in analysis: unknown case!" << endl;
			exit(1);
			break;
			// something seriously wrong here!
			}

		count_multisites++;
		}	/*end loop over sites*/

	set(_bialsites, locus, count_bialsites);
	set(_multisites, locus, count_multisites);
	if (outgroup)
		set(_sfout, locus, count_sfA + count_sfB);
	else
		set(_sfout, locus, count_sfout);
	set(_sfA, locus, count_sfA);
	set(_sfB, locus, count_sfB);
	if (outgroup)
		{
		set(_sxA, locus, count_sxA + count_sxAfB);
		set(_sxB, locus, count_sxB + count_sxBfA);
		}
	else
		{
		set(_sxA, locus, count_sxA);
		set(_sxB, locus, count_sxB);
		set(_sxAfB, locus, count_sxAfB);
		set(_sxBfA, locus, count_sxBfA);
		}
	set(_ss, locus, count_ss);
	set(_Wald, locus, (double) WaldWolfowitz(Walds, 0));
	}	/*end of compute_polyl*/

/**********************************************************************************/

/*From Miguel Navascues*/
// calculation of p-values has been removed

// Function to calculate number of runs or its standardize (i.e. normalized) value
// ARGUMENTS:
//   number_of_elements: length of the sequence
//   IOseq: pointer to a set of integers that represents the sequence to test
//          this sequence should contain exclusively "1" and "0" (ones and zeros)
//	 std: which statistic should be returned: 0=number of runs; 1=standardize number of runs
//        with any other value the number of runs is returned
//        (note that the number of runs is an integer but will be returned as a double)
double WaldWolfowitz(const wald_type & IOseq, int stats)
	{
	int runs = 1;
	int n0 = 0, n1 = 0;

	for (size_t i=0; i<IOseq.size(); i++)
		{
		IOseq[i] ? n1++ : n0++;
		
		if (i>0 && IOseq[i]!=IOseq[i-1]) 
			runs++;
		}

	if (n0<5 || n1<5) 
		return -9;
	
	if (n0+n1 != int(IOseq.size())) 
		cout << "n0+n1<>n \n";
	
	if (runs > int(IOseq.size())) 
		cout << "runs>n \n";
	
	if (stats==1)
		{
		const double mean = 1.0 + 
			( (double(n0) * double(n1) * 2.0) / double(IOseq.size()));
		const double variance = 
			( (mean-1) * (mean-2) ) / (double(IOseq.size()-1));
		return ( double(runs) - mean ) / sqrt(variance);
		}		
	else
		return runs;
	}



double pearson_corr_pi(const SingleStats & statsA, 
	const SingleStats & statsB) 
	{
	const int nloc = statsA[_pi].size();

    double ay=0.0, ax=0.0;

    for (int j=0; j<nloc; j++)
   		{
        ax += statsA[_pi][j];
        ay += statsB[_pi][j];
    	}

    ax /= nloc;
    ay /= nloc;

    double syy=0.0, sxy=0.0, sxx=0.0;

    for (int j=0; j<nloc; j++)
		{
        const double xt = statsA[_pi][j]-ax;
        const double yt = statsB[_pi][j]-ay;
        sxx += xt*xt;
        syy += yt*yt;
        sxy += xt*yt;
		}
    return sxy/sqrt(sxx*syy);
	}




