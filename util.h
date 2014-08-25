#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <limits>


template<class T>
class CachedValue
	{
protected:
	bool _ready;
	T _value;

public:
	CachedValue()
		: _ready(false)
		{}

	void reset()
		{
		_ready = false;
		}

	bool ready() const
		{
		return _ready;
		}

	void set_ready()
		{
		_ready = true;
		}

	const T & operator=(const T & value)
		{
		set_ready();
		_value = value;
		return _value;
		}

	const T & operator()(const T & value)
		{
		set_ready();
		_value = value;
		return _value;
		}

	operator const T&() const
		{
		return _value;
		}

	T& get() 
		{
		return _value;
		}
	};


template<class ITER>
struct pair_iter : public std::pair<ITER, ITER>
	{
	typedef typename ITER::value_type single_v_type;
	typedef std::pair<single_v_type *, single_v_type *> reference;

	pair_iter(const ITER & i1, const ITER & i2)
		: std::pair<ITER, ITER>(i1, i2)
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
