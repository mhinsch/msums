#include <iostream>
#include <vector>
#include <utility>

#include "alleleset.h"
#include "stats_multi.h"

using namespace std;

void add(AlleleSet<char> & set, int n0, int n1)
	{
	for (int i=0; i<n0; i++)
		set.add(0);

	for (int i=0; i<n1; i++)
		set.add(1);
	}

int main()
	{
	vector<AlleleSet<char> > p1, p2;

	p1.resize(60); p2.resize(60);

	add(p1[2], 0, 3);
	add(p1[11], 1, 2);
	add(p1[17], 0, 3);
	add(p1[24], 2, 1);
	add(p1[26], 1, 2);
	add(p1[31], 2, 1);
	add(p1[34], 1, 2);
	add(p1[35], 2, 1);
	add(p1[39], 1, 2);
	add(p1[44], 2, 1);
	add(p1[50], 2, 1);
	add(p1[55], 0, 3);
	add(p1[59], 0, 3);

	add(p2[0], 2, 1);
	add(p2[5], 0, 3);
	add(p2[10], 0, 3);
	add(p2[19], 1, 2);
	add(p2[24], 1, 2);
	add(p2[31], 0, 3);
	add(p2[35], 2, 1);
	add(p2[48], 1, 2);
	add(p2[52], 0, 3);
	add(p2[56], 1, 2);

	for (size_t i=0; i<p1.size(); i++)
		{
		if (p1[i].size() == 0)
			add(p1[i], 3, 0);
		p1[i].do_count();

		if (p1[i].size() == 0)
			{
			cerr << "size 0 at " << i << endl;
			exit(1);
			}

		if (p1[i].count(0) + p1[i].count(1) != 3)
			{
			cerr << "wrong number at " << i << endl;
			exit(1);
			}
		}
	
	for (size_t i=0; i<p2.size(); i++)
		{
		if (p2[i].size() == 0)
			add(p2[i], 3, 0);
		p2[i].do_count();

		if (p2[i].size() == 0)
			{
			cerr << "size 0 at " << i << endl;
			exit(1);
			}

		if (p2[i].count(0) + p2[i].count(1) != 3)
			{
			cerr << "wrong number at " << i << endl;
			exit(1);
			}
		}
	
	pair<double, double> res = W_navascues_01(
		p1.begin(), p1.end(),
		p2.begin(), p2.end());

	cout << res.first << " " << res.second << "\n";

	PolyNavascues poly;
	poly.compute(
		p1.begin(), p1.end(),
		p2.begin(), p2.end());

	cout << poly.count_rf << " " << poly.count_rs << "\n";
	}
