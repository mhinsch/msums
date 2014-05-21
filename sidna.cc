#include <cstdlib>

#include <string>
#include <vector>
#include <iostream>

#include <boost/dynamic_bitset.hpp>

#include "spmain.h"
#include "sidna.h"

using namespace std;

/** Defines the structure used to store results from shared sites, 
  	fixed sites per locus. */
struct biallelic_sites
	{
	int n;
	char b1;
	char b2;
	};


/** Checks whether number of alleles at @site is > 2. */
bool check_multi(const vector<string> & sequences, 
	int site, biallelic_sites & result)
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


int pair_dif(const vector<string> & data, int nspoly)
	{
	int diff = 0;
	for(size_t h1=0; h1<(data.size()-1); h1++) 
		for(size_t  h2=h1+1; h2<data.size(); h2++) 
			for(int s=0; s<nspoly; s++) 
				if(data[h1][s] != data[h2][s]) 
					diff++;

	return diff;
	}

int avg_pair_div(const vector<string> & seqAlhs, const vector<string> & seqBlhs,
	int nspolyl)
	{
	int count_dif = 0;
	for(size_t h1=0; h1<seqAlhs.size(); h1++) 
		for(size_t h2=0; h2<seqBlhs.size(); h2++) 
			for(int s=0; s<nspolyl; s++) 
				if(seqAlhs[h1][s] != seqBlhs[h2][s]) count_dif++;

	return count_dif;
	}

int segr(const vector<string> & data, int nspolyl)
	{
	/*now computes theta for species A */
	int count_segr = 0;
	for(int s=0; s<nspolyl; s++) 
		{	/*starts loop over sites at locus l*/
		const char first_base = data[0][s];
		for(size_t h=1; h<data.size(); h++) 
			{	/*starts loop over haplotypes at species A*/
			if(data[h][s] != first_base) 
				{	/*current site is segregating*/
				count_segr++;
				break;
				}
			}	/*end loop over haplotypes in A*/
		}	/*end loop over sites*/
	
	return count_segr;
	}


typedef vector<char> wald_type;
double WaldWolfowitz(const wald_type & IOseq, int stats);
float DTajima(int nseq, int nsegsites, float avgpairdif);



void SingleStats::compute_pi_theta(const size_t locus, 
	const vector<string> & seqlhs, int nspolyl)
	{
	const int nseq = seqlhs.size();
	const int nSites = seqlhs[0].size();

	/*compute pi for species A */
	const float sumpairdif = pair_dif(seqlhs, nspolyl) / float(nSites);
	const int count_pair = seqlhs.size() * (seqlhs.size() - 1) / 2;
	set(_sumpairdif, locus, sumpairdif);
	set(_pi, locus, sumpairdif / count_pair);
	
	// theta
	const int count_segr = segr(seqlhs, nspolyl);
	set(_theta, locus, count_segr / _watt(nseq) / float(nSites));
	set(_D, locus, DTajima(nseq, count_segr, get(_pi, locus)*nSites));	/*Warning Pi is multiplied by number of sites*/
	}


void PairStats::compute_div_fst(const size_t locus,
	const vector<string> & seqAlhs, const vector<string> & seqBlhs,
	const SingleStats & statsA, const SingleStats & statsB,
	int nspolyl)
	{
	const int nseqA = seqAlhs.size();
	const int nseqB = seqBlhs.size();
	const int nSites = seqAlhs[0].size();

	/*now computes average pairwise divergence between species A and B*/
	const int count_pairAB = nseqA * nseqB;
	const float pair_div = avg_pair_div(seqAlhs, seqBlhs, nspolyl) / float(nSites);

	set(_d, locus, pair_div / float(count_pairAB));
	set(_dn, locus, 
		get(_d, locus) - (statsA.get(_pi, locus) + statsB.get(_pi, locus))/2.0);

	const int count_pair_piT = 
		nseqA * (nseqA - 1) / 2 + nseqB * (nseqB - 1) / 2 + count_pairAB;
	
	const float piT = 
		(statsA.get(_sumpairdif, locus) + statsB.get(_sumpairdif, locus) + 
		pair_div) / (float) count_pair_piT;

	if (piT >= 1.0e-7) 
		set(_FST, locus, 
			(piT - (statsA.get(_pi, locus)+statsB.get(_pi, locus))/2.0)/piT);
	}

void PairStats::compute_polyl(const size_t locus,
	const vector<string> & seqAlhs, const vector<string> & seqBlhs, 
	int nspolyl, const string & seqOls, bool outgroup)
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
	set(_Wald, locus, (float) WaldWolfowitz(Walds, 0));
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


float DTajima(int nseq, int nsegsites, float avgpairdif)
	{
	float a1=0, a2=0;
		
	for(int i=1; i<nseq; i++) 
		{
		a1 += (float) 1.0 / i; 
		a2 += (float) 1.0 / i / i;
		}

	const float b1 = (float) (nseq+1.0) / 3.0 / (nseq-1.0);
	const float b2 = (float) 2.0 * (nseq*nseq+nseq+3.0) / 9.0 / nseq / (nseq-1.0);
	const float c1 = b1 - (1.0/a1);
	const float c2 = b2 - ((nseq+2.0)/a1/nseq) + (a2/a1/a1);
	const float e1 = c1 / a1;
	const float e2 = c2 / (a1*a1+a2);

	const double Dtemp = nsegsites == 0 ? 
		0.0F : /*Warning put a value of 0 if no segregating sites!!!!!!*/
		(double) (avgpairdif-((double) nsegsites/(double)a1))
			/ sqrt(e1*nsegsites+e2*nsegsites*(nsegsites-1)); 
	
	return (float) Dtemp;
	}	/*end of procedure DTajima*/        


float pearson_corr_pi(const SingleStats & statsA, 
	const SingleStats & statsB) 
	{
	const int nloc = statsA[_pi].size();

    float ay=0.0, ax=0.0;

    for (int j=0; j<nloc; j++)
   		{
        ax += statsA[_pi][j];
        ay += statsB[_pi][j];
    	}

    ax /= nloc;
    ay /= nloc;

    float syy=0.0, sxy=0.0, sxx=0.0;

    for (int j=0; j<nloc; j++)
		{
        const float xt = statsA[_pi][j]-ax;
        const float yt = statsB[_pi][j]-ay;
        sxx += xt*xt;
        syy += yt*yt;
        sxy += xt*yt;
		}
    return sxy/sqrt(sxx*syy);
	}




