#ifndef MSUMSOPTIONS_H
#define MSUMSOPTIONS_H

#include <string>

#include "sputil.h"

#include "options.h"

class MSUMSOptions
	{
public:
	string inputfilename;
	string statfilename;
	bool printPerLocus;

	typedef OpList::Collector OLIter;

	OpList popList, statList, groupStatList;

private:
	bool printHelp;
	bool listStats;
	
	OLIter keepS, dropS, keepP, dropP, groupSt;

	OptionParser p;
	UnaryOption<string> o_ifn;
	UnaryOption<string> o_sfn;
	SwitchOption o_ppl;
	SwitchOption o_hlp;
	SwitchOption o_list;
	NaryOption<string, OLIter> o_ds;
	NaryOption<string, OLIter> o_ks;
	NaryOption<string, OLIter> o_dp;
	NaryOption<string, OLIter> o_kp;
	NaryOption<string, OLIter> o_ms;
	
public:
	MSUMSOptions()
		: inputfilename("spinput.txt"), statfilename("ABCstat.txt"),
		printPerLocus(false), printHelp(false), listStats(false),
		keepS(statList, "+"), dropS(statList, "-"), 
		keepP(popList, "+"), dropP(popList, "-"),
		groupSt(groupStatList, ""), 
		o_ifn(p, 'i', "init", inputfilename), o_sfn(p, 'o', "output", statfilename),
		o_ppl(p, 'l', "print_per_locus", printPerLocus),
		o_hlp(p, 'h', "help", printHelp),
		o_list(p, 'a', "list_stats", listStats),
		o_ds(p, 's', "dropStats", dropS), o_ks(p, 'S', "keepStats", keepS),
		o_dp(p, 'p', "dropPops", dropP), o_kp(p, 'P', "keepPops", keepP),
		o_ms(p, 'm', "multiStats", groupSt)
		{}

	void do_print_help();
	void do_list_stats();

	void process(int argc, char * argv[])
		{
		try {
		p.parse(argv+1, argc-1);
		} catch (ArgException & e)
			{
			error(string("Error in main: unknown command line argument: ") +
				e.what());
			}

		if (printHelp)
			{
			do_print_help();
			exit(0);
			}

		if (listStats)
			{
			do_list_stats();
			exit(0);
			}
		}
	};


#endif	// MSUMSOPTIONS_H
