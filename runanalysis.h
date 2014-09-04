#ifndef RUNANALYSIS_H
#define RUNANALYSIS_H

#include <iostream>
#include <vector>
#include <string>

#include "anafunctors.h"
#include "sample.h"
#include "pairsample.h"
#include "setup.h"

using namespace std;

void analyse(ostream & abcstatfile, int dataset,
	vector<vector<StrSample> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SingleStats*> & single_stats,
	const vector<PairStats*> & pair_stats,
	const vector<GroupStats*> & group_stats,
	const vector<vector<size_t> > & groups);

void analyse_aggr(ostream & abcstatfile, int dataset,
	vector<vector<StrSample> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SingleStats*> & single_stats,
	const vector<PairStats*> & pair_stats,
	const vector<GroupStats*> & group_stats,
	const vector<vector<size_t> > & groups);

#endif	// RUNANALYSIS_H
