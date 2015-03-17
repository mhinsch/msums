#ifndef OPTIONS_H
#define OPTIONS_H

#include <vector>
#include <string>

using namespace std;

class ArgException : public exception
	{
protected:
	string _msg;
public:
	ArgException(const string & msg)
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

	virtual ~ArgException() throw()
		{}
	};


class OptionParser
	{
public:
	class Option
		{
	public:
		Option(OptionParser & p)
			{
			p.add(*this);
			}

		virtual char ** match(char ** argv, int argc) = 0;
		};

	void add(Option & option)
		{
		_options.push_back(&option);
		}

	void parse(char ** argv, int argc)
		{
		char ** end = argv + argc;
		for (char ** c=argv; c<end; c++)
			{
			size_t p;
			for (p=0; p<_options.size(); p++)
				{
				char ** r = _options[p]->match(c, end-c);
				if (r)
					{
					c = r-1;
					break;
					}
				}

			if (p==_options.size())
				throw ArgException(*c);
			}
		}

protected:

	vector<Option *> _options;
	};


class SwitchOption : public OptionParser::Option
	{
public:
	SwitchOption(OptionParser & p, const char flag, const string & name, bool & var,
		bool to=true)
		: OptionParser::Option(p), _var(var), _to(to), _name(name), _flag(flag)
		{}

	char ** match(char ** argv, int argc)
		{
		if (argc<1)
			return 0;

		const string arg(argv[0]);
		if (arg == string("-") + _flag ||
			arg == "--" + _name)
			{
			_var = _to;
			return argv + 1;
			}
		
		// no match
		return 0;
		}

protected:
	bool & _var, _to;
	string _name;
	char _flag;
	};

template<class T>
class UnaryOption : public OptionParser::Option
	{
public:
	UnaryOption(OptionParser & p, const char flag, const string & name, T & var)
		: OptionParser::Option(p), _var(var), _name(name), _flag(flag)
		{}

	char ** match(char ** argv, int argc)
		{
		if (argc<2)
			return 0;

		const string arg(argv[0]);
		if (arg == string("-") + _flag ||
			arg == "--" + _name)
			{
			_var = lexical_cast<T>(argv[1]);
			return argv + 2;
			}
		
		// no match
		return 0;
		}

protected:
	T & _var;
	string _name;
	char _flag;
	};

template<typename T, typename I>
class NaryOption : public OptionParser::Option
	{
public:
	NaryOption(OptionParser & p, const char flag, const string & name, 
		I & iter)
		: OptionParser::Option(p), _iter(iter), _name(name), _flag(flag)
		{}

	char ** match(char ** argv, int argc)
		{
		if (argc<2)
			return 0;

		const string arg(argv[0]);
		if (arg == string("-") + _flag ||
			arg == "--" + _name)
			{
			int i=1;

			_iter.begin();
			
			for (; i<argc; i++)
				{
				if (argv[i][0] == '-')
					break;

				_iter(lexical_cast<T>(argv[i]));
				}			

			return argv + i;
			}
		
		// no match
		return 0;
		}

protected:
	I & _iter;
	string _name;
	char _flag;
	};

// Argument with n parameters that can re-occur several times in the cmdl.
// Single occurences are added as a new line to a matrix preceded by an operator 
// that is set in the Collector object. Enables e.g. '--add a b c --remove x y z' 
// to be stored in a single data structure.
class OpList
	{
public:
	struct Collector
		{
		vector<vector<string> > & _list;
		string _op;

		// empty op strings will not be added to the list
		Collector(OpList & list, const string & op = "")
			: _list(list.getList()), _op(op)
			{}

		void begin()
			{
			_list.push_back(vector<string>());
			if (_op != "")
				_list.back().push_back(_op);
			}

		void operator()(const string & str)
			{
			//cerr << "ol: " << str << endl;
			_list.back().push_back(str);
			}
		};

	vector<vector<string> > & getList()
		{
		return _list;
		};

protected:
	vector<vector<string> > _list;
	};


#endif	// OPTIONS_H
