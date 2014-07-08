
#include <iterator>

template<class SET_ITER>
int sum_pair_diff(
	const SET_ITER & start1, const SET_ITER & stop1, 
	const SET_ITER & start2, const SET_ITER & stop2
	)
	{
	typedef std::iterator_traits<SET_ITER>::value_type::const_iterator all_iter;
	int count = 0;
	// iterate through sites on both allele count lists
	for(SET_ITER set1=start1, set2=start2; 
		set1!=stop1 && set2!=stop2; 
		set1++, set2++)
		// compare allele counts
		for(all_iter s1=set1->begin(); s1!=set1->end(); s1++) 
		for(all_iter s2=set2->begin(); s2!=set2->end(); s2++) 
			if (s1->first != s2->first)
				count += s1->second * s2->second;

	return count;
	}


template<class SET_ITER, class SEQ_ITER>
void compute_polyl(
	const SET_ITER & start1, const SET_ITER & stop1, 
	const SET_ITER & start2, const SET_ITER & stop2, 
	const SEQ_ITER & start_out, const SEQ_ITER & stop_out)
	{
	const int nSites = s1.size();

	int count_bialsites = 0;
	int count_multisites = 0;
	int count_sfA = 0;
	int count_sfB = 0;
	int count_sfout = 0;
	int count_sxA = 0;
	int count_sxB = 0;
	int count_sxAfB = 0;
	int count_sxBfA = 0;
	int count_ss = 0;

	typedef std::iterator_traits<SEQ_ITER>::value_type state_type;
	/*vector of runs of fixed and polymorphic sites according to Miguel Navascues*/
	std::vector<state_type> Walds();

	SEQ_ITER s_out = start_out;

	for (SET_ITER s1=start1, s2=start2; 
		s1!=stop1 && s2!=stop2; 
		s1++, s2++, s_out++) 
		{	/*starts loop over polymorphic sites at locus l*/
		int n_all1 = s1->size();
		int n_all2 = s2->size();

		if (n_all1 > 2 || n_all2 > 2) 
			{	/*multiple hits*/
			count_multisites++;
			continue;	/*skip to next site*/
			}

		// only fixed or biallelic arrives here

		const state_type outgroup = *s_out;

		// now test whether multiple hits and if not 
		// determine the type of polymorphism
		const int code = n_all1 + n_all2*2;

		switch (code)
			{
		case 3:	// a==1, b==1
			// all identical
			if (*s1[0]==*s2[0] && *s1[0]==outgroup)
				continue;
			// all different => multiallelic
			if (*s1[0]!=*s2[0] && *s1[0]!=outgroup && *s2[0]!=outgroup)
				break;

			// two alleles
			count_bialsites++;

			if (*s1[0] == *s2[0])
				/*site fixed along outgroup*/
				count_sfout++;
			else	// exactly one of a,b identical to outgroup
				{	
				Walds[s] = true;
				*s1[0] == outgroup ?
					count_sfB++ :	// base in B is derived
					count_sfA++; 	// base in A is derived
				}

			continue;
		case 4: // a==2, b==1
			if ( (*s2[0]==*s1[0] || *s2[0]==*s1[1]) &&
				(outgroup==*s2[0] || outgroup==*s1[0] || outgroup==*s1[1]) ) 
				{	/* fixed in B for one of the two bases in A*/
				Walds[s] = false;
				count_bialsites++;

				if (outgroup == *s2[0]) 	/*unique derived polymorphism in A*/
					count_sxA++;
				else					/*unique ancestral polymorphism in A*/
					count_sxAfB++;

				continue;
				} 
			break;
		case 5: // a==1, b==2
			// fixed in A for one out of B and one of them equal to outgroup
			if ( (*s1[0]==*s2[0] || *s1[0]==*s2[1]) &&
				(outgroup==*s1[0] || outgroup==*s2[0] || outgroup==*s2[1]) )
				{
				Walds[s] = false;
				count_bialsites++;

				if (outgroup == *s1[0])	// unique derived polymorphism in B
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
			if (  ( (*s1[0]==*s2[0] && *s1[1]==*s2[1]) ||
					(*s1[0]==*s2[1] && *s1[1]==*s2[0]) )  &&
				(outgroup==*s1[0] || outgroup==*s1[1])  )
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
	set(_Wald, locus, (double) WaldWolfowitz(Walds.begin(), Walds.end(), 0));
	}	/*end of compute_polyl*/



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
