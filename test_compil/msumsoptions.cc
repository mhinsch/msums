#include <iostream>


#include "msumsoptions.h"

using namespace std;


void MSUMSOptions::do_print_help()
	{
	cout <<
	"analms [-i FILE] [-o FILE] [-h] [-l] [[-s|-S STATLIST]...] "
	"[[-p|-P POPLIST]...]\n\n";
	cout <<
	"-i, --init FILE       Use this initialization file. [spinput.txt]\n"
	"-o, --output FILE     Use this output file. [ABCstat.txt]\n"
	"-h, --help            Print the help text and exit.\n"
	"-a, --list_stats      Print list of available stats and exit.\n"
	"-l, --print_per_locus Print stats per locus, don't print aggregate values.\n"
	"                      [false]\n"
	"-S, --keepStats LIST  Print stats in LIST. LIST is a space-separated list of\n"
	"                      stat names. 'all' selects all stats.\n"
	"-s, --dropStats LIST  Do *not* print stats in LIST (see above).\n"
	"-P, --keepPops LIST   Show populations or pairs in LIST. Populations can be\n"
	"                      specified as single number (e.g. '2'), or range \n"
	"                      (e.g. '1-5'). Pairs are specifies as p1xp2 (e.g. '1x2').\n"
	"                      'all' selects all populations, 'allxall' all pairs.\n"
	"-p, --dropPops LIST   Do *not* show populations/pairs in LIST (see above).\n"
	"-m, --multiStats LIST Show multi-population statistics in LIST. Use e.g. as:\n"
	"                      -m fst 1x2x3 3x4x5 -m f4 1x4x5x6\n";
	}

void MSUMSOptions::do_list_stats()
	{
	cout << "sorry, not implemented yet :-(\n";
	}
