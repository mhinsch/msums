#include <cstdlib>

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "setup.h"
#include "statheader.h"
#include "sputil.h"
#include "msumsoptions.h"
#include "runanalysis.h"
#include "configfile.h"
#include "msdatafile.h"

using namespace std;


typedef GroupStatHandler<vector<StrSample*> > gs_handler_t;

void select_stats(const vector<vector<string> > & sel, 
	const vector<string> & all, set<string> & use);
void select_pops(const vector<vector<string> > & sel,
	vector<bool> & pops, vector<vector<bool> > & pairs);
void create_group_stats(
	const vector<vector<string> > & sel, gs_handler_t & gsh);



int main(int argc,char *argv[])
	{
	MSUMSOptions opt;

	opt.process(argc, argv);


// **** read config file

	ifstream inp(opt.inputfilename.c_str());
	if (!inp.good())
		error("Error in main: cannot open config file " + opt.inputfilename);

	ConfigFile conf;
	try {
	conf.read(inp);
	} catch (exception & e) {error(e.what());}
	
	conf.dump(cout);



// **** now we know #pops we can process command line arguments

	// ** pops

	vector<bool> show_pops(conf.n_pops(), true);
	vector<vector<bool> > show_pairs;
	for (size_t p=0; p<conf.n_pops(); p++)
		show_pairs.push_back(vector<bool>(conf.n_pops(), true));
	// process command line selection of pops/pairs
	select_pops(opt.popList.getList(), show_pops, show_pairs);

	// ** stats

	SSHandler ss_handler;
	PSHandler ps_handler;

	// names of *all* stats (single and pair)
	vector<string> all_stat_names(
		ss_handler.names().begin(), ss_handler.names().end());
	all_stat_names.insert(all_stat_names.end(),
		ps_handler.names().begin(), ps_handler.names().end());

	set<string> selected_stats;
	// process command line selection of stats
	select_stats(opt.statList.getList(), all_stat_names, selected_stats);
	
	// set the selected ones to active
	ss_handler.activate(selected_stats);
	ps_handler.activate(selected_stats);

	const vector<SingleStats*> & single_stats = ss_handler.active();
	const vector<PairStats*> & pair_stats = ps_handler.active();	

	gs_handler_t gs_handler;
	create_group_stats(opt.groupStatList.getList(), gs_handler);
	const vector<gs_handler_t::analysis_t*> & group_stats = gs_handler.stats();
	const vector<gs_handler_t::group_t> & groups = gs_handler.groups();
	

// **** prepare header of output file

	ofstream abcstatfile(opt.statfilename.c_str());
	if (!abcstatfile.good())
		error("Error in main: cannot open output file " + opt.statfilename);

	try {
	write_ABCstat_header(abcstatfile, 
		show_pops, show_pairs, 
		ss_handler.active_names(), ps_handler.active_names(), 
		gs_handler.names(), gs_handler.groups(),
		!opt.printPerLocus);
	} catch (exception & e)
		{error(e.what());}


// **** prepare data containers

	/** Sequence by species x locus x sample. */
	vector<vector<StrSample> > sequences;

	size_t maxNSeq = 0;

	// building array for sequence data
	sequences.resize(conf.n_pops(), vector<StrSample>(0));
	for (size_t p=0; p<sequences.size(); p++)
		{
		sequences[p].resize(conf.n_loci(), StrSample());
		for (size_t l=0; l<sequences[p].size(); l++)
			{
			sequences[p][l].resize(conf.n_sequences()[p][l], Sequence());
			maxNSeq = max(maxNSeq, conf.n_sequences()[p][l]);
			}
		}

	/** Outgroup data. */
	vector<Sequence> outgroup(conf.n_loci(), Sequence());
	
	ifstream input(conf.datafilename().c_str());
	if (!input.good())
		error("Error in main: cannot open data file " +
			conf.datafilename());


// **** start reading datasets and calculating
// **** dataset = replicate, i.e. 1 dataset = loci x pops x samples

	for(size_t dataset=0; dataset<conf.n_datasets(); ++dataset) 
		{ 
		const int imain = conf.n_datasets()/100;
		if (imain>0 && (dataset % imain)==0) 
			{
			string t_str;
			get_time_short(t_str);
			cerr << "\nReplicate: " << dataset << " Time: " << t_str;
			}


		for (size_t p=0; p<sequences.size(); p++)
			for (size_t l=0; l<sequences[p].size(); l++)
				sequences[p][l].reset();

		try {
		read_dataset(input, dataset, 
			conf.n_sequences(), conf.n_sites(), sequences, outgroup);
		} catch (exception & e)
			{error(e.what());}

		// new data set
            // the function statheader ends with 'endl', so no need to introduce a newline for the first dataset
		if (dataset > 0)
            {
			abcstatfile << "\n";
            }

		if (opt.printPerLocus)
			analyse(abcstatfile, dataset, sequences, show_pops, show_pairs,
				single_stats, pair_stats, group_stats, groups);
		else
			analyse_aggr(abcstatfile, dataset, sequences, show_pops, show_pairs,
				single_stats, pair_stats, group_stats, groups);
		}

	return 0;
	}    /*end Main*/


void select_stats(
	const vector<vector<string> > & sel, 
	const vector<string> & all, set<string> & use)
	{
	for (size_t s=0; s<sel.size(); s++)
		{
		cerr << "sS: " << sel[s][0];
		const bool to = sel[s][0][0] == '+';
		for (size_t t=1; t<sel[s].size(); t++)
			{
			cerr << ", " << sel[s][t];
			if (sel[s][t] == "all")
				{
				for (size_t s=0; s<all.size(); s++)
					if (to)
						use.insert(all[s]);
					else
						use.erase(all[s]);
				continue;
				}

			if (to)
				use.insert(sel[s][t]);
			else
				use.erase(sel[s][t]);
			}
		cerr << endl;
		}
	}

void create_group_stats(
	const vector<vector<string> > & sel, gs_handler_t & gsh)
	{
	vector<size_t> selected;
	// go through all selection sequences
	for (size_t s=0; s<sel.size(); s++)
		{
		const vector<string> & list = sel[s];
		const string & stat_name = list[0];

		// go through list of pops
		for (size_t t=1; t<list.size(); t++)
			{
			selected.clear();
			const string & str = list[t];
			back_insert_iterator<vector<size_t> > sel_inserter(selected);

			try {
			splitStr(str, 'x', sel_inserter);
			} catch(exception & e)
				{error(string("Error in parameter: '") + str + 
					"' is not a valid group");}
			
			gsh.add(stat_name, selected);
			}
		}
	}

void checkPopIndex(size_t i, size_t n)
	{
	if (i >= n)
		error(string("Error in parameter: population index ") + 
			lexical_cast<string>(i) + " out of range");
	}

void select_pops(const vector<vector<string> > & sel,
	vector<bool> & pops, vector<vector<bool> > & pairs)
	{
	// go through all selection sequences
	for (size_t s=0; s<sel.size(); s++)
		{
		// prefix determines whether this adds or removes pops/pairs
		const bool to = sel[s][0][0] == '+';
		
		// go through list of operations (skipping prefix)
		for (size_t t=1; t<sel[s].size(); t++)
			{
//			cerr << "-pP " << sel[s][t] << endl;
			const string & str = sel[s][t];
			// all pops
			if (str == "all")
				{
				for (size_t p=0; p<pops.size(); p++)
					pops[p] = to;
				continue;
				}
			// all pairs
			if (str == "allxall")
				{
				for (size_t p=0; p<pairs.size()-1; p++)
					for (size_t p2=p+1; p2<pairs[p].size(); p2++)
						pairs[p][p2] = to;
				continue;
				}

			size_t c;
			for (c=0; c<str.size(); c++)
				{
				if (str[c] == '-')	// pop range
					{
					size_t p1, p2;

					try {
					splitStr(str, c, p1, p2);
					} catch(exception & e)
						{error(string("Error in parameter: '") + str + 
							"' is not a valid population range");}

					checkPopIndex(p1, pops.size());
					checkPopIndex(p2, pops.size());

					for (size_t p=p1; p<=p2; p++)
						{
						pops[p] = to;
						for (size_t pp=0; pp<pops.size(); pp++)
							pairs[min(p, pp)][max(p, pp)] = pops[p] && pops[pp];
						}
					break;
					}
				if (str[c] == 'x') // pair
					{
					int p1, p2;

					try {
					splitStr(str, c, p1, p2);
					} catch(exception & e)
						{error(string("Error in parameter: '") + str + 
							"' is not a valid pair");}

					checkPopIndex(p1, pops.size());
					checkPopIndex(p2, pops.size());

					pairs[min(p1,p2)][max(p1,p2)] = to;
					break;
					}
				}

			if (c == str.size()) // no '-' or 'x' found => single pop
				{
				size_t p;
				try {
				p = lexical_cast<size_t>(str);
				} catch(exception & e)
					{error(string("Error in parameter: '") + str + 
						"' is not a valid population or pair");}

				checkPopIndex(p, pops.size());

				pops[p] = to;
				for (size_t pp=0; pp<pops.size(); pp++)
					pairs[min(p, pp)][max(p, pp)] = pops[p] && pops[pp];
				}
			}
		}
	}

