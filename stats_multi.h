#ifndef STATS_MULTI_H
#define STATS_MULTI_H

#include <cmath>
#include <cstdlib>

#include <iterator>
#include <utility>
#include <vector>

#include "patterson.h"


template<class SET_ITER>
int sum_pair_diff(
	const SET_ITER & start1, const SET_ITER & stop1, 
	const SET_ITER & start2, const SET_ITER & stop2
	)
	{
	int count = 0;
	// iterate through sites on both allele count lists
	for(SET_ITER set1=start1, set2=start2; 
		set1!=stop1 && set2!=stop2; 
		set1++, set2++)
		// compare allele counts
		for(size_t s1=0; s1<set1->size(); s1++) 
		for(size_t s2=0; s2<set2->size(); s2++) 
			if (s1 != s2)
				count += (*set1)[s1] * (*set2)[s2];

	return count;
	}


struct PolyL
	{
	int bialsites, multisites, sfout, sfA;
	int sfB, sxA, sxB, sxAfB, sxBfA;
	int ss;
	double wald;

	template<class SET_ITER, class SEQ_ITER>
	void compute(
		const SET_ITER & start1, const SET_ITER & stop1, 
		const SET_ITER & start2, const SET_ITER & stop2, 
		const SEQ_ITER & start_out, const SEQ_ITER & stop_out)
		{
		typedef typename std::iterator_traits<SEQ_ITER>::value_type state_type;
		/*vector of runs of fixed and polymorphic sites according to Miguel Navascues*/
		std::vector<char> Walds;

		SEQ_ITER s_out = start_out;

		// go through all sites
		for (SET_ITER s1=start1, s2=start2; 
			s1!=stop1 && s2!=stop2; 
			s1++, s2++, s_out++) 
			{
			const int n_all1 = s1->n_alleles();
			const int n_all2 = s2->n_alleles();

			if (n_all1 > 2 || n_all2 > 2) 
				{	/*multiple hits*/
				multisites++;
				continue;	/*skip to next site*/
				}

		// only fixed or biallelic arrives here

			const state_type outgroup = *s_out;

			// now test whether multiple hits and if not 
			// determine the type of polymorphism
			const int code = n_all1 + n_all2*2;

			// should be checked, but we know they have to be valid
			state_type s1_0 = s1->find();
			state_type s2_0 = s2->find();
			// these might return invalid, but switch over code will take 
			// care of that
			state_type s1_1 = s1->find(s1_0+1);
			state_type s2_1 = s2->find(s2_0+1);

			switch (code)
				{
			case 3:	// a==1, b==1
				// all identical
				if (s1_0==s2_0 && s1_0==outgroup)
					continue;
				// all different => multiallelic
				if (s1_0!=s2_0 && s1_0!=outgroup && s2_0!=outgroup)
					break;
				// from here on 2 alleles
				
				bialsites++;

				if (s1_0 == s2_0)
					// both derived
					sfout++;
				else	// exactly one of a,b identical to outgroup
					{	
					Walds[s] = 1;
					s1_0 == outgroup ?
						sfB++ :	// base in B is derived
						sfA++; 	// base in A is derived
					}

				continue;
			case 4: // a==2, b==1
				if ( (s2_0==s1_0 || s2_0==s1_1) &&
					(outgroup==s2_0 || outgroup==s1_0 || outgroup==s1_1) ) 
					{	/* fixed in B for one of the two bases in A*/
					Walds[s] = 0;
					bialsites++;

					if (outgroup == s2_0) 	/*unique derived polymorphism in A*/
						sxA++;
					else					/*unique ancestral polymorphism in A*/
						sxAfB++;

					continue;
					} 
				break;
			case 5: // a==1, b==2
				// fixed in A for one out of B and one of them equal to outgroup
				if ( (s1_0==s2_0 || s1_0==s2_1) &&
					(outgroup==s1_0 || outgroup==s2_0 || outgroup==s2_1) )
					{
					//Walds[s] = false;
					bialsites++;

					if (outgroup == s1_0)	// unique derived polymorphism in B
						sxB++;
					else 					// unique ancestral polymorphism in B
						sxBfA++;

					continue;
					} 
				// base fixed in A is different from both in B
				// or outgroup different from all
				break;
			case 6:	// a==2, b==2
				// shared polymorphism and outgroup present
				if (  ( (s1_0==s2_0 && s1_1==s2_1) ||
						(s1_0==s2_1 && s1_1==s2_0) )  &&
					(outgroup==s1_0 || outgroup==s1_1)  )
					{
					ss++;
					//Walds[s] = false;
					bialsites++;
					continue;
					} 
				/* polymorphism for different pairs of bases in A and B*/
				break;
			default:
				std::cerr << "Error in analysis: unknown case!" << std::endl;
				exit(1);
				break;
				// something seriously wrong here!
				}

			multisites++;
			}	/*end loop over sites*/

	/*	if (outgroup)
			{
			sfout = sfA + sfB;
			sxA += sxAfB;
			sxB += sxBfA;
			}
	*/
		//wald = (double) WaldWolfowitz(Walds.begin(), Walds.end(), 0);
		}	/*end of compute_polyl*/
	};

// calculate W statistics following Navascues et al. 2014
// ONLY VALID FOR MS-STYLE DATA!
// (states 0 - ancestral, 1 - derived)
template<class SET_ITER>
std::pair<double, double> navascues_W_01(
	const SET_ITER & start1, const SET_ITER & stop1, 
	const SET_ITER & start2, const SET_ITER & stop2) 
	{
	static std::vector<int> runs1, runs2;

	runs1.clear();
	runs2.clear();

	int n_poly = 0;
	int runl1 = 0, runl2 = 0;

	for (SET_ITER s1=start1, s2=start2; 
		s1!=stop1 && s2!=stop2; 
		s1++, s2++) 
		{
		const int n1 = s1->n_alleles(), n2 = s2->n_alleles();
		if (n1==1 && n2==1)
			continue;

		bool c;
		
		if (n1==2 && n2==2)
			{
			const int na_1 = s1->count(0), nd_1 = s1->count(1);
			const int na_2 = s2->count(0), nd_2 = s2->count(1);
			const double f1 = nd_1/double(na_1+nd_1);
			const double f2 = nd_2/double(na_2+nd_2);
			if (f1==f2) continue;

			c = f1>f2;
			}
		else 	// one pop fixed, one poly
			c = n1 == 1;

		n_poly++;

		if (c)	// 1 for pop1, 0 for pop2
			{
			runl2++;
			runs1.push_back(runl1);
			runl1 = 0;
			}
		else	// 0 for pop1, 1 for pop2
			{
			runl1++;
			runs2.push_back(runl2);
			runl2 = 0;
			}
		}

	runs1.push_back(runl1);
	runs2.push_back(runl2);

	double kmin2 = n_poly - 2;

	double bynp1 = 1.0/runs1.size();
	double wx2s1 = 0;
	for (size_t i=0; i<runs1.size(); i++)
		wx2s1 += fabs(runs1[i]/kmin2 - bynp1);
	wx2s1 /= 2;
	

	double bynp2 = 1.0/runs2.size();
	double wx1s2 = 0;
	for (size_t i=0; i<runs2.size(); i++)
		wx1s2 += fabs(runs2[i]/kmin2 - bynp2);
	wx1s2 /= 2;

	return std::pair<double, double>(wx2s1, wx1s2);
	}

// classify sites following Navascues et al. 2014
// F == 0; S == 3
template<class SET>
int classify(const SET & s1, const SET & s2)
	{
	const bool p1 = s1.n_alleles()>1;
	const bool p2 = s2.n_alleles()>1;

	if (p1 && p2)
		return 3;

	if (!p1 && !p2)
		return s1.product(s2) > 0 ? 1 : 0;

	return 2;
	}


template<class SET_ITER>
std::pair<int, int> navascues_R(
	const SET_ITER & start1, const SET_ITER & stop1, 
	const SET_ITER & start2, const SET_ITER & stop2) 
	{
	// F == 0b00 == 0u; S == 0b11 == 3u
	SET_ITER s1=start1, s2=start2; 

	int c;

	// throw away leading 1s
	for (; s1!=stop1 && s2!=stop2 && (c=classify(*s1, *s2))==1; 
		s1++, s2++)
		{}

	if (s1==stop1 || s2==stop2) return std::pair<int, int>(0, 0);

	// found first polymorphic locus
	
	bool rf = c == 0;
	bool rs = c == 3;
	int count_rf = 1;
	int count_rs = 1;

	for (s1++, s2++; s1!=stop1 && s2!=stop2; 
		s1++, s2++) 
		{
		c = classify(*s1, *s2);
		if (c==1) continue;

		std::cout << (c==0);
		if ((c == 0) != rf)	// rf state changes => new run
			{
			count_rf++;
			rf = ! rf;
			}

		if ((c == 3) != rs)	// rs state changes => new run
			{
			count_rs++;
			rs = ! rs;
			}
		}

	return std::pair<int, int>(count_rf, count_rs);
	}


/*From Miguel Navascues*/
// calculation of p-values has been removed

// Function to calculate number of runs or its standardize (i.e. normalized) value
// @stats which statistic should be returned: 0=number of runs; 1=standardize number of runs
//        with any other value the number of runs is returned
//        (note that the number of runs is an integer but will be returned as a double)
// **CURRENTLY ASSUMES 0/1 DATA!**
template<class SEQ_ITER>
double WaldWolfowitz(const SEQ_ITER & start, const SEQ_ITER & stop, int stats)
	{
	int runs = 1;
	int n0 = 0, n1 = 0;
	SEQ_ITER prev = start;

	int size = 0;

	for (SEQ_ITER s=start; s!=stop; s++)
		{
		*s ? n1++ : n0++;
		
		if (s!=start && *s!=*prev) 
			runs++;

		prev = s;
		size++;
		}

	if (n0<5 || n1<5) 
		return -9;
	
	if (stats==1)
		{
		const double mean = 1.0 + 
			( (double(n0) * double(n1) * 2.0) / double(size));
		const double variance = 
			( (mean-1) * (mean-2) ) / (double(size-1));
		return ( double(runs) - mean ) / sqrt(variance);
		}		
	else
		return runs;
	}

template<class SET_ITER>
double patterson_D(
	const SET_ITER & start1, const SET_ITER & stop1,
	const SET_ITER & start2, const SET_ITER & stop2)
	{
	int n = 0;
	double d = 0;

	for (SET_ITER s1=start1, s2=start2; s1!=stop1; s1++, s2++)
		{
		n++;
		d += D_est(s1->count(0), s1->count(1), s2->count(0), s2->count(1));
		}

	return d / n;
	}


template<class SET_ITER>
double patterson_f3
	(const SET_ITER & start1, const SET_ITER & stop1,
	 const SET_ITER & start2, const SET_ITER & stop2,
	 const SET_ITER & start3, const SET_ITER & stop3)
	{
	int n = 0;
	double f = 0;

	for (SET_ITER s1=start1, s2=start2, s3=start3;
		s1!=stop1; s1++, s2++, s3++)
		{
		f += f3_est(
			s1->count(0), s1->count(1),
			s2->count(0), s2->count(1),
			s3->count(0), s3->count(1));
		n++;
		}

	return f / n;
	}

template<class SET_ITER>
double patterson_f4
	(const SET_ITER & start1, const SET_ITER & stop1,
	 const SET_ITER & start2, const SET_ITER & stop2,
	 const SET_ITER & start3, const SET_ITER & stop3,
	 const SET_ITER & start4, const SET_ITER & stop4)
	{
	int n = 0;
	double f = 0;

	for (SET_ITER s1=start1, s2=start2, s3=start3, s4=start4;
		s1!=stop1; s1++, s2++, s3++, s4++)
		{
		f += f4_est(
			s1->count(0), s1->count(1),
			s2->count(0), s2->count(1),
			s3->count(0), s3->count(1),
			s4->count(0), s4->count(1));
		n++;
		}

	return f / n;
	}

#endif	// STATS_MULTI_H
