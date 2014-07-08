#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <limits>

template<class T, class U>
class Harmonic_Rec
	{
public:
	typedef T data_t;
	
	T operator()()
		{
		return T(0);
		}

	T operator()(T last, U curr)
		{
		return last + T(1)/T(curr);
		}
	};

template<class T, class U>
class Harmonic_2_Rec
	{
public:
	typedef T data_t;
	
	T operator()()
		{
		return T(0);
		}

	T operator()(T last, U curr)
		{
		return last + T(1)/T(curr*curr);
		}
	};

template<class OP>
class RecursiveSeries
	{
public:
	typedef typename OP::data_t data_t;
	
protected:
	std::vector<data_t> _data;
	OP _op;

public:
	RecursiveSeries()
		: _data(1)
		{
		_data[0] = _op();
		}
	
	template<class I>
	data_t operator()(I n)
		{
		for (size_t i=_data.size(); i<=n; i++)
			_data.push_back(_op(_data.back(), i));
		
		return _data[n];
		}
	};

typedef RecursiveSeries<Harmonic_Rec<double, int> > Harmonic;
typedef RecursiveSeries<Harmonic_2_Rec<double, int> > Harmonic_2;

inline double undefined()
	{
	return std::numeric_limits<double>::quiet_NaN();
	}

template<class I>
I num_pairs(I n) { return n*(n-1)/2; }

class InvalidEnsure : public std::exception
	{
public:
	const char * what()
		{
		return "No operations on const object possible!";
		}
	};

template<class CLASS>
struct EnsureCalled
	{
	template<void (CLASS::*FUNC)()>
	static bool was_called()
		{
		static bool called = false;

		if (called) return true;

		return !(called = true);
		}

	template<void (CLASS::*FUNC)()>
	bool did(const CLASS & object) const
		{
		if (was_called<FUNC>())	// all is well
			return false;

		// can't call function on const object!
		throw InvalidEnsure();
		}

	template<void (CLASS::*FUNC)()>
	bool did(CLASS & object)() const
		{
		if (was_called<FUNC>())
			return false;

		object.*FUNC();
		return true;
		}
	};

template<class CLASS>
struct WrapCached
	{
	struct Base
		{
		bool called;
		const CLASS & object;

		Base(bool c, CLASS & o)
			: called(c), object(o)
			{}
		};

	template<class T, T (CLASS::* FUNC)()>
	struct CachedValue : public Base
		{
		T value;

		Cached(CLASS & obj)
			: Base(false, obj)
			{}

		T operator()()
			{
			if (called)
				return value;
			called = true;
			return value = object.*FUNC();
			}
		};
	};


template<class ITER>
struct pair_iter : public std::pair<ITER, ITER>
	{
	typedef ITER::value_type single_v_type;
	typedef std::pair<single_v_type *, single_v_type *> reference;

	pair_iter(const ITER & i1, const ITER & i2)
		: pair(i1, i2)
		{}

	pair_iter & operator++()
		{
		this->first++;
		this->second++;
		return *this;
		}

	pair_iter operator++(int)
		{
		pair_iter ret(*this);
		++*this;
		return ret;
		}
	
	reference operator*()
		{
		return reference(&(*this->first), &(*this->second));
		}
	};


#endif	//UTIL_H
