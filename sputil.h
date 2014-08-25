#ifndef SPUTIL_H
#define SPUTIL_H


#include <ctime>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <string>
#include <exception>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


#define STRFY(s) #s
#define TOSTR(s) STRFY(s)
#define ERR_LOC __FILE__ ":" TOSTR(__LINE__)

using boost::lexical_cast;

/*************************************************************************/

inline double pow_2(double x)
	{
	return (double) x*(double)x;
	}

/*************************************************************************/

inline void get_time_short(std::string & smess)
	/*get the current time and put only the day number and time in smess*/
	{
	smess = lexical_cast<std::string>(time(NULL));
	}


inline void error(const std::string & msg)
	{
	std::cerr << "\n" << msg << std::endl;
	exit(1);
	}

template<typename T>
void splitStr(const std::string & str, int pos, T & x1, T & x2)
	{
	x1 = lexical_cast<T>(str.substr(0, pos));
	x2 = lexical_cast<T>(str.substr(pos+1, str.size()-pos-1));
	}

template<typename ITER>
void splitStr(const std::string & str, char sep, ITER & iter)
	{
	typedef typename ITER::container_type::value_type value_type;
	size_t last = 0;
	for (size_t i=0; i<str.size(); i++)
		if (str[i] == sep || i==str.size()-1)
			{
			*iter = lexical_cast<value_type>(str.substr(last, i-last));
			last = i + 1;
			}
	}

inline void skip_space(std::istream & inp, std::string & str)
	{
	while(getline(inp, str) && boost::all(str, boost::is_space()));
	}

template<typename T> 
bool get_value (std::istream & inp_file, T & value)
	{
	static std::string str;

	skip_space(inp_file, str);
	if (!inp_file) return false;

	std::istringstream tmp_s;
	tmp_s.str(str);
	tmp_s >> value;
	if (tmp_s.fail()) return false;

	return true;
	}


class SPIOException : public std::exception
	{
protected:
	std::string _msg;
public:
	SPIOException(const std::string & msg)
		: _msg(msg)
		{}

	const char * what() const throw()
		{
		return _msg.c_str();
		}

	void print(std::ostream & out) const
		{
		out << _msg;
		}

	virtual ~SPIOException() throw()
		{}
	};



#endif // SPUTIL_H
