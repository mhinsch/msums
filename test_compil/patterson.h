#ifndef PATTERSON_H
#define PATTERSON_H

// giving s is redundant but faster
inline double h_est(int n1, int n2, int s)
	{
	return n1 * n2 / double(s * (s-1));
	}

inline double freq(int n1, int n2)
	{
	return n1 / double(n1 + n2);
	}

inline double f2_est(int n1_a, int n2_a, int n1_b, int n2_b)
	{
	const int sa = n1_a + n2_a, sb = n1_b + n2_b;
	const double a = n1_a / double(sa);
	const double b = n1_b / double(sb);

	return (a-b)*(a-b) - 
		h_est(n1_a, n2_a, sa)/sa - h_est(n1_b, n2_b, sb)/sb;
	}

inline double f3_est(int n1_a, int n2_a, int n1_b, int n2_b, int n1_c, int n2_c)
	{
	const int sc = n1_c + n2_c;
	const double a = freq(n1_a, n2_a);
	const double b = freq(n1_b, n2_b);
	const double c = n1_c / double(sc);

	return (c-a)*(c-b) - h_est(n1_c, n2_c, sc)/sc;
	}

inline double f4_est(int n1_a, int n2_a, int n1_b, int n2_b, 
	int n1_c, int n2_c, int n1_d, int n2_d)
	{
	return (freq(n1_a, n2_a) - freq(n1_b, n2_b)) *
		(freq(n1_c, n2_c) - freq(n1_d, n2_d));
	}

inline double D_est(int n1_a, int n2_a, int n1_b, int n2_b)
	{
	const int sa = n1_a + n2_a, sb = n1_b + n2_b;
	const double a = n1_a / double(sa);
	const double b = n1_b / double(sb);
	const double h_a = h_est(n1_a, n2_a, sa);
	const double h_b = h_est(n1_b, n2_b, sb);

	return (a-b)*(a-b) + h_a*(1 - 1.0/sa) + h_b*(1 - 1.0/sb);
	}

#endif	// PATTERSON_H
