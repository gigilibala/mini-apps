#include <iostream>
#include <list>
#include <vector>

class Entry {
public:
	const char* name;
	double samples;
	int count;
//	std::list<double> samples;

	Entry(const char* name) : name(name), samples(0.0) , count(0) { };
	~Entry() { };
};

/* Used for timing benchmarks */
class Benchmark {

private:
	std::vector<Entry*> entries;

public:
	Benchmark() { };
	~Benchmark() { };

	/* Adds an entry for timing and returns the id to that entry. */
	int    add_entry(const char* name);
	void   add_sample(int entry_id, double sample);

	double get_sum(int entry_id);
	double get_mean(int entry_id);
	
	double get_low_mean_percentile(int percentile);
};
