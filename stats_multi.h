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

	PolyL()
		: bialsites(0), multisites(0), sfout(0), sfA(0), sfB(0), sxA(0), sxB(0),
		  sxAfB(0), sxBfA(0), ss(0)
	{}

	template<class SET_ITER>
	void compute(
		const SET_ITER & start1, const SET_ITER & stop1, 
		const SET_ITER & start2, const SET_ITER & stop2)
		{
		typedef typename std::iterator_traits<SET_ITER>::value_type::state_t state_type;

		const bool has_outgroup = true;

		// go through all sites
		for (SET_ITER s1=start1, s2=start2; 
			s1!=stop1 && s2!=stop2; 
			s1++, s2++) 
			{
			const int n_all1 = s1->n_alleles();
			const int n_all2 = s2->n_alleles();
			// cout << "NA: " << n_all1 << " " << n_all2 << "\n";

			if (n_all1 > 2 || n_all2 > 2) 
				{	/*multiple hits*/
				multisites++;
				continue;	/*skip to next site*/
				}

		// only fixed or biallelic arrives here

			// hard-coded for now
			const state_type outgroup = 0;

			// now test whether multiple hits and if not 
			// determine the type of polymorphism
			const int code = n_all1 + n_all2*2;

			// should be checked, but we know they have to be valid
			const state_type sA_0 = s1->find();
			const state_type sB_0 = s2->find();
			// these might return invalid, but switch over code will take 
			// care of that
			const state_type sA_1 = s1->find(sA_0+1);
			const state_type sB_1 = s2->find(sB_0+1);
			// cout << "s: " << code << " (" << n_all1 << " " << n_all2 << "); " << sA_0 << " " << sA_1 << ", " << sB_0 << " " << sB_1 << ", " << outgroup << "\n";

			switch (code)
				{
			case 3:	// a==1, b==1
				// all identical
				if (sA_0==sB_0 && sA_0==outgroup)
					continue;
				// cout << ".";
				// all different => multiallelic
				if (sA_0!=sB_0 && sA_0!=outgroup && sB_0!=outgroup)
					break;
				// from here on 2 alleles
			
				// cout << "~";	
				bialsites++;

				if (sA_0 == sB_0)
					// both derived
					sfout++;
				else	// exactly one of a,b identical to outgroup
					sA_0 == outgroup ?
						sfB++ :	// base in B is derived
						sfA++; 	// base in A is derived

				continue;
				// cout << "!!";	
			case 4: // a==2, b==1
				if ( (sB_0==sA_0 || sB_0==sA_1) &&
					(outgroup==sB_0 || outgroup==sA_0 || outgroup==sA_1) ) 
					{	/* fixed in B for one of the two bases in A*/
					bialsites++;

					// cout << "+";

					if (outgroup == sB_0) 	/*unique derived polymorphism in A*/
						sxA++;
					else					/*unique ancestral polymorphism in A*/
						sxAfB++;

					continue;
					} 
				// cout << "-";
				break;
			case 5: // a==1, b==2
				// fixed in A for one out of B and one of them equal to outgroup
				if ( (sA_0==sB_0 || sA_0==sB_1) &&
					(outgroup==sA_0 || outgroup==sB_0 || outgroup==sB_1) )
					{
					bialsites++;

					if (outgroup == sA_0)	// unique derived polymorphism in B
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
				if (  ( (sA_0==sB_0 && sA_1==sB_1) ||
						(sA_0==sB_1 && sA_1==sB_0) )  &&
					(outgroup==sA_0 || outgroup==sA_1)  )
					{
					ss++;
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

		// cout << "S" << ss << " " << multisites << "\n";
		if (has_outgroup)
			{
			sfout = sfA + sfB;
			sxA += sxAfB;
			sxB += sxBfA;
			}
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

	if (n_poly <= 2)
		return std::pair<double, double>(0, 0);

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

		//std::cout << (c==0);
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
