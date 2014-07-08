#include "stats.h"


Harmonic harmonic;
Harmonic_2 harmonic_2;
	

double DTajima(int n_sequences, int n_segsites, double theta_pi)
	{
	if (n_segsites==0)
		return 0;

	const double a1 = harmonic(n_sequences-1);
	const double a2 = harmonic_2(n_sequences-1);
	const double b1 = double(n_sequences+1) / (3 * (n_sequences-1));
	const double b2 = 
		2 * double(n_sequences*n_sequences + n_sequences + 3) 
			/ 
		(9 * n_sequences * (n_sequences-1));
	const double c1 = b1 - (1/a1);
	const double c2 = b2 - (n_sequences+2)/(a1*n_sequences) + a2/(a1*a1);
	const double e1 = c1 / a1;
	const double e2 = c2 / (a1*a1+a2);

	return (theta_pi - n_segsites/a1)
				/ 
			sqrt(e1*n_segsites + e2*n_segsites*(n_segsites-1)); 
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
		sqrt(uF*n_mutations+ vF*n_mutations*n_mutations);
	return FStar;
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
