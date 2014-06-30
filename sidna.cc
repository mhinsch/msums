#include <cstdlib>

#include <vector>
#include <iostream>

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

Harmonic harmonic();
Harmonic_2 harmonic_2();
	

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

int count_singletons(const Sample & sequences, const size_t polym_sites)
	{
	bialellic_sites bs;
	int count = 0;

	for (size_t site=0; site<polym_sites; site++)
		if (check_multi(sequences, site, bs))
			count++;

	return count;
	}

int pair_dif(const Sample & data, int n_poly_sites)
	{
	int diff = 0;
	for(size_t h1=0; h1<(data.size()-1); h1++) 
		for(size_t  h2=h1+1; h2<data.size(); h2++) 
			for(int s=0; s<n_poly_sites; s++) 
				if(data[h1][s] != data[h2][s]) 
					diff++;

	return diff;
	}

int avg_pair_div(const Sample & sample1, const Sample & sample2,
	int n_poly_sites)
	{
	int count_dif = 0;
	for(size_t h1=0; h1<sample1.size(); h1++) 
		for(size_t h2=0; h2<sample2.size(); h2++) 
			for(int s=0; s<n_poly_sites; s++) 
				if(sample1[h1][s] != sample2[h2][s]) count_dif++;

	return count_dif;
	}

// number of polymorphic sites in sequence
int segr(const Sample & sample, int n_poly_sites)
	{
	int count_segr = 0;
	for(int site=0; site<n_poly_sites; site++) 
		{
		const char first_base = sample[0][site];
		for(size_t h=1; h<sample.size(); h++) 
			{
			if(sample[h][site] != first_base) 
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


float DTajima(int n_sequences, int n_segsites, float avgpairdif)
	{
	const float a1 = harmonic(n_sequences-1);
	const float a2 = harmonic_2(n_sequences-1);
	const float b1 = (n_sequences+1.0) / 3.0 / (n_sequences-1.0);
	const float b2 = 2.0 * (n_sequences*n_sequences+n_sequences+3.0) / 
		9.0 / n_sequences / (n_sequences-1.0);
	const float c1 = b1 - (1.0/a1);
	const float c2 = b2 - ((n_sequences+2.0)/a1/n_sequences) + (a2/a1/a1);
	const float e1 = c1 / a1;
	const float e2 = c2 / (a1*a1+a2);

	const double Dtemp = nsegsites == 0 ? 
		0.0F : /*Warning put a value of 0 if no segregating sites!!!!!!*/
		(double) (avgpairdif-((double) n_segsites/(double)a1))
			/ sqrt(e1*n_segsites+e2*n_segsites*(n_segsites-1)); 
	
	return (float) Dtemp;
	}	/*end of procedure DTajima*/        


void SingleStats::compute_pi_theta(const size_t locus, 
	const Sample & sample, int n_polym_sites)
	{
	const int n_sequences = sample.size();
	const int n_sites = sample[0].size();

	/*compute pi for species A */
	const float sum_pair_dif = pair_dif(sample, n_polym_sites) / float(n_sites);
	const int count_pair = n_sequences * (n_sequences - 1) / 2;
	set(_sumpairdif, locus, sum_pair_dif);
	set(_pi, locus, sum_pair_dif / count_pair);
	
	// theta
	const int count_segr = segr(sample, n_polym_sites);
	set(_theta, locus, count_segr / harmonic(n_sequences-1) / float(n_sites));
	set(_D, locus, DTajima(n_sequences, count_segr, get(_pi, locus)*n_sites));	/*Warning Pi is multiplied by number of sites*/
	}


double c_sub_n(int n_sequences)
	{
	if (n_sequences <= 2)
		return 1.0;

	const double a = harmonic(n_sequences - 1);
	const double c = 2 * (n_sequences * a - 2 * (n_sequences - 1));
	return c / ((n_sequences - 1) * (n_sequences - 2));
	}

double d_sub_n(int n_sequences)
	{
	const double a_n_plus1 = harmonic(n_sequences);
	const int nsm1 = n_sequences - 1;
	const int nsm2 = n_sequences - 2;
	double d = c_sub_n(n_sequences) + double(nsm2) / (nsm1 * nsm1);
	d += (2.0/nsm1) * (1.5 - (2*a_n_plus1 - 3)/nsm2 -  1.0/n_sequences);
	return d;
	}


double fu_li_Dstar(int n_sequences, int n_singletons, int n_mutations)
	{
	const double a = harmonic(n_sequences-1);
	const double b = harmonic_2(n_sequences-1);
	const double d = d_sub_n(n_sequences);

	const int nsm1 = n_sequences - 1;
	const double nsbynsm1 = double(n_sequences) / nsm1;
	const double vD = 
		(nsbynsm1*nsbynsm1*b + a*a*d - 2*(n_sequences * a * (a + 1))/(nsm1 * nsm1))
			/
		(a*a + b);

	const double uD = nsbynsm1*(a - nsbynsm1) - vD;

	const double DStar = 
		(nsbynsm1*n_mutations - a*n_singletons)
		/
		sqrt(uD*n_mutations + vD*n_mutations*n_mutations);
	return DStar;
	}

double fu_li_Fstar(int n_sequences, int n_singletons, int n_mutations,
	double pi)
	{
	const double a = harmonic(n_sequences - 1);
	const double a_n_plus1 = harmonic(n_sequences);
	const double b = harmonic_2(n_sequences - 1);

	const int ns_2 = n_sequences * n_sequences;
	const int ns_3 = ns_2 * n_sequences;

	// according to libsequence (v 1.8):
	//vF is taken from the correction published by
	//Simonsen et al.  (1995) Genetics 141: 413, eqn A5
	const double vF = 
		(
		(2*ns_3 + 110*ns_2 - 255*n_sequences + 153) / (9 * ns_2 * (n_sequences - 1))
		+ 
		2*(n_sequences - 1)*a/ns_2 - 8*b/n_sequences
		)
		/ 
		(a*a + b);

	const double uF = 
		(4*ns_2 + 19*n_sequences + 3 - 12*(n_sequences + 1)*a_n_plus1)
		/ 
		(3 * n_sequences * (n_sequences - 1) * a)
		- vF;

	const double FStar = 
		(pi - double(n_sequences - 1)/n_sequences*n_singletons)
		/
		sqrt(uF*NumMut + vF*n_mutations*n_mutations);
	return FStar;
	}

void SingleStats::compute_fu_li(const size_t locus,
	const Sample & sample, int n_polym_sites)
	{
	const double pi = get(_pi, locus);
	const int n_mutations = segr(sample, n_polym_sites);
	const int n_singletons = count_singletons(sample, n_polym_sites);
	
	const double d_star = fu_li_Dstar(sample.size(), n_singletons, n_mutations);
	const double f_star = fu_li_Fstar(sample.size(), n_singletons, n_mutations, pi);

	set(_Dstar, locus, d_star);
	set(_Fstar, locus, f_star);
	}

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




