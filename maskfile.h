#ifndef MASKFILE_H
#define MASKFILE_H

#include "sputil.h"

void read_mask(istream & inp, vector<vector<bool> > & mask)
	{
	string str;
	
	while(true)
		{
		skip_space(inp, str);
		
		if (! inp.good())
			break;

		mask.push_back(vector<bool>());
		mask.back().resize(str.size());

		for (size_t i=0; i<str.size(); i++)
			mask.back()[i] = str[i] == '1';
		}
	}

#endif	// MASKFILE_H
