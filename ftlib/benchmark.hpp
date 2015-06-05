#include <iostream>
#include <list>


using namespace std;


class Benchmark{

	list<double> samples;

	Benchmark();
	~Benchmark();

	void add_sample(double sample);

	double get_low_mean_percentile(int percentile);
};
