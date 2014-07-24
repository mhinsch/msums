#ifndef RUNANALYSIS_H
#define RUNANALYSIS_H

#include <iostream>
#include <vector>
#include <string>

#include "setup.h"

using namespace std;

void analyse(ostream & abcstatfile, int dataset,
	vector<vector<StrSample> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SSHandler::analysis_t*> & single_stats,
	const vector<PSHandler::analysis_t*> & pair_stats,
	const vector<gs_handler_t::analysis_t*> & group_stats,
	const vector<vector<size_t> > & groups);

void analyse_aggr(ostream & abcstatfile, int dataset,
	vector<vector<StrSample> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SSHandler::analysis_t*> & single_stats,
	const vector<PSHandler::analysis_t*> & pair_stats,
	const vector<gs_handler_t::analysis_t*> & group_stats,
	const vector<vector<size_t> > & groups);

#endif	// RUNANALYSIS_H
