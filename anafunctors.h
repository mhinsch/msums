#ifndef ANAFUNCTORS_H
#define ANAFUNCTORS_H


template<class SAMPLE>
struct AnalysisBase
	{
	virtual double analyse(const SAMPLE & s) const = 0;
	double operator()(const SAMPLE & s) const = 0
		{ return this->analyse(s); }
	};

template<class SAMPLE, class STAT>
class Analysis : public AnalysisBase
	{
protected:
	STAT SAMPLE::* _stat_p;

public:
	Analysis(STAT SAMPLE::* stat_p)
		: _stat_p(stat_p)
		{}

	double analyse(const SAMPLE & s) const
		{
		return (s.*stat_p)();
		}
	};

template<class SAMPLE, class STAT>
Analysis<SAMPLE, STAT> * create_analysis(STAT SAMPLE::* stat_p)
	{
	return new Analysis<SAMPLE, STAT>(stat_p);
	}


struct Aggregate
	{
protected:
	double _mean, _sqr_sum, _std;
	size_t _n;

public:
	Aggregate()
		: _n(0), _mean(0), _sqr_sum(0), _std(0)
		{}

	void operator()(double value)
		{
		if (value != value)
			return;

		_n++;
		_mean += value;
		_sqr_sum += value * value;
		}

	void analyse()
		{
		if (n<2)
			return;
		
		const double ssq = _sqr_sum - _mean*_mean/n;
		_std = sqrt(ssq / n-1);

		_mean /= n;
		}

	double mean() const
		{ return _mean; }
	double std() const
		{ return _std; }
	};


#endif	// ANAFUNCTORS_H
