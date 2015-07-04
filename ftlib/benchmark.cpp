#include "benchmark.hpp"

int Benchmark::add_entry(const char* name) {
	Entry* entry = new Entry(name);
	entries.push_back(entry);
	return entries.size()-1;
}

void Benchmark::add_sample(int entry_id, double sample) {
	entries[entry_id]->samples += sample;
	entries[entry_id]->count ++;
}

double Benchmark::get_sum(int entry_id) {
	return entries[entry_id]->samples;
}

double Benchmark::get_mean(int entry_id) {
	return (get_sum(entry_id) / entries[entry_id]->count);
}
