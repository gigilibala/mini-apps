#include "benchmark.hpp"
#include <mpi.h>


int Benchmark::add_entry(const char* name) {
	Entry* entry = new Entry(name);
	entries.push_back(entry);
	return entries.size()-1;
}

void Benchmark::add_sample(int entry_id, double sample) {
	entries[entry_id]->samples += sample;
	entries[entry_id]->count ++;
}

void Benchmark::start_timing(int entry_id) {
	entries[entry_id]->t1 = MPI_Wtime();
}

void Benchmark::end_timing(int entry_id) {
	add_sample(entry_id, MPI_Wtime() - entries[entry_id]->t1);
}

double Benchmark::get_sum(int entry_id) {
	return entries[entry_id]->samples;
}

double Benchmark::get_mean(int entry_id) {
	return (get_sum(entry_id) / entries[entry_id]->count);
}
