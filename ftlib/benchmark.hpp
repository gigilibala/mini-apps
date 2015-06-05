#include <iostream>

using namespace std;


class Benchmark{

	list<double> samples;

	Benchmark();
	~Benchmark();

	void add_sample(double sample);

	double get_mean_percentile(int percentile);
};
