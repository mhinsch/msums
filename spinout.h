#ifndef SPINOUT_H
#define SPINOUT_H


#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>

#include "sidna.h"
#include "spmain.h"


using namespace std;

class SPIOException : public exception
	{
protected:
	string _msg;
public:
	SPIOException(const string & msg)
		: _msg(msg)
		{}

	const char * what() const throw()
		{
		return _msg.c_str();
		}

	void print(ostream & out) const
		{
		out << _msg;
		}

	virtual ~SPIOException() throw()
		{}
	};


/** Reads the first non-whitespace string. */
void skip_space(istream & inp, string & str);

void check(istream & inp, ostream & err, const string & str);

template<typename T> 
bool get_value (istream & inp_file, T & value)
	{
	static std::string str;

	skip_space(inp_file, str);
	if (!inp_file) return false;

	istringstream tmp_s;
	tmp_s.str(str);
	tmp_s >> value;
	if (tmp_s.fail()) return false;

	return true;
	}

void get_initial_conditions
	(const string & inputfile, vector<vector<int> > & nSeq, vector<int> & nSites, 
	 int & ndatasets, string & datafile);

void print_initial_conditions(ostream & out, 
	const string & inputfile, vector<vector<int> > nSeq, vector<int> nSites, 
	long ndatasets, const string & datafile);

/** Reads datasetfile and extracts DNA sequences from all species and 
	the outgroup.
	@param inp file to read data from.
	@param dataset number of current data set.
	@param nSeq #sequences[population][locus].
	@param nSites number of sites[locus].
	@param seqhs sequences[population][locus][number].
	@param seqOls outgroup data[locus].
	@param nspolyl 
	*/
void get_dataset(istream & inp, 
	long dataset, const vector<vector<int> > & nSeq, const vector<int> & nSites, 
	vector<vector<vector<string> > > & seqhs, vector<string> & seqOls, 
	vector<int> & nspolyl);

/** Initialize the output file where writing statistics of shared and fixed 
	polymorphism to outputfile one line for each replicate dataset. */
void initialize_write_ABCstat(ostream & file, 
	const vector<bool> & pops, const vector<vector<bool> > & pairs, 
	const vector<bool> sstats, const vector<bool> pstats,
	const StatDict & dict, bool writeAggr);


template<typename STATLIST>
void writeResults(ostream & out, 
	const STATLIST & results, const vector<bool> & stats)
	{
	for (size_t i=0; i<stats.size(); i++)
		{
		if (!stats[i])
			continue;
		out << '\t' << results[typename STATLIST::idx_type(i)].mean() <<
			'\t' << results[typename STATLIST::idx_type(i)].std();
		}
	}

template<typename STATLIST>
void writeLResults(ostream & out,
	const STATLIST & results, const vector<bool> & stats, int locus)
	{
	for (size_t i=0; i<stats.size(); i++)
		{
		if (!stats[i])
			continue;
		out << '\t' << results[typename STATLIST::idx_type(i)][locus];
		}
	}


#endif	//SPINOUT_H
