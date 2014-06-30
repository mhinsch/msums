#ifndef SPUTIL_H
#define SPUTIL_H


#include <ctime>
#include <cstdlib>

#include <iostream>

#include <boost/lexical_cast.hpp>


#define STRFY(s) #s
#define TOSTR(s) STRFY(s)
#define ERR_LOC __FILE__ ":" TOSTR(__LINE__)

using boost::lexical_cast;
using namespace std;

/*************************************************************************/

inline double pow_2(double x)
	{
	return (double) x*(double)x;
	}

/*************************************************************************/

inline void get_time_short(string & smess)
	/*get the current time and put only the day number and time in smess*/
	{
	smess = lexical_cast<string>(time(NULL));
	}


inline void error(const string & msg)
	{
	cerr << "\n" << msg << endl;
	exit(1);
	}


#endif // SPUTIL_H
