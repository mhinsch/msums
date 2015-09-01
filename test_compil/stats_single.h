#ifndef STATS_H
#define STATS_H

#include <algorithm>
#include <cmath>

#include <iostream>
using namespace std;

#include "mathutil.h"

typedef RecursiveSeries<Harmonic_Rec<double, int> > Harmonic;
typedef RecursiveSeries<Harmonic_2_Rec<double, int> > Harmonic_2;
	
extern Harmonic harmonic;
extern Harmonic_2 harmonic_2;


// number of polymorphic sites
template<class SET_ITER>
int n_segregating_sites(const SET_ITER & start, const SET_ITER & stop)
	{
	int n = 0;
	for (SET_ITER set = start; set != stop; set++)
		if (set->n_alleles() > 1)
			n++;

	return n;
	}


// sum of pairwise differences over all sites
template<class SET_ITER>
int pairwise_differences(const SET_ITER & start, const SET_ITER & stop)
	{
	int sum = 0;

	for (SET_ITER set=start; set!=stop; set++)
		sum += set->pairwise_difference();
	
	return sum;
	}


// number of singletons over all sites
template<class SET_ITER>
int count_singletons(const SET_ITER & start, const SET_ITER & stop)
	{
	int count = 0;

	for (SET_ITER set=start; set!=stop; set++)
		count += set->n_singletons();

	return count;
	}


template<class SET_ITER, class SEQ_ITER>
int count_singletons_for_seq(
	const SEQ_ITER & seq_start, const SEQ_ITER & seq_stop,
	 // !! implicitly assumes same size for seq and sets
	const SET_ITER & set_start)
	{
	int count = 0;

	SET_ITER s = set_start;
	// go through sites
	for (SEQ_ITER i=seq_start; i!=seq_stop; i++, s++)
		if (s->is_singleton(*i))
			count++;

	return count;
	}

inline double theta_pi(int n_sequences, int sum_pair_diff)
	{
	return sum_pair_diff * 2.0 / (n_sequences * (n_sequences-1));
	}


inline double theta_W(int n_sequences, int n_segr)
	{
	return n_segr / harmonic(n_sequences-1);
	}

// Tajima's D
double DTajima(int n_sequences, int n_segsites, double theta_pi);


inline double c_sub_n(int n_sequences)
	{
	if (n_sequences <= 2)
		return 1.0;

	const double a = harmonic(n_sequences - 1);
	const double c = 2 * (n_sequences * a - 2 * (n_sequences - 1));
	return c / ((n_sequences - 1) * (n_sequences - 2));
	}


inline double d_sub_n(int n_sequences)
	{
	const double a_n_plus1 = harmonic(n_sequences);
	const int nsm1 = n_sequences - 1;
	const int nsm2 = n_sequences - 2;
	double d = c_sub_n(n_sequences) + double(nsm2) / (nsm1 * nsm1);
	d += (2.0/nsm1) * (1.5 - (2*a_n_plus1 - 3)/nsm2 -  1.0/n_sequences);
	return d;
	}


double fu_li_Fstar(int n_sequences, int n_singletons, int n_mutations, double pi);


double fu_li_Dstar(int n_sequences, int n_singletons, int n_mutations);


template<class ITER>
double R2(const ITER & start, const ITER & stop, double pi, int S)
	{
    if (S == 0 || start==stop) return undefined();

    double sm = 0.0;
	int c = 0;

    for (ITER i=start; i!=stop; i++, c++)
            sm += (*i - pi/2.0) * (*i - pi/2.0);

    sm = sqrt(sm/c) / S;

    if (sm < 1.0E-15)
            sm = 0.0;

    return sm;
	}

#endif	// STATS_H
