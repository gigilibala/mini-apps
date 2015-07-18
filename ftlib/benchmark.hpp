/*
 * Copyright (c) 2013-2014 University of Alabama at Birmingham
 *                         Department of Computer and Information Sciences
 *                         All rights reserved.
 * $COPYRIGHT$
 *
 * Author: Amin Hassani (ahassani@cis.uab.edu)
 *******************************************************************************
 * 
 * 
 */

#include <iostream>
#include <list>
#include <vector>
#include <mpi.h>
#include <stdio.h>

#ifndef FTLIB_BENCHMARK_HPP
#define FTLIB_BENCHMARK_HPP

/* Used for timing benchmarks */
class BenchmarkEntry {
public:
	const char* name;
	double samples;
	int count;
	double t1;
//	std::list<double> samples;
	char out_str[100];
	BenchmarkEntry(const char* name) : name(name), samples(0.0) , count(0) { };
	~BenchmarkEntry() { };

	void   start_timing() { t1 = MPI_Wtime(); };
	void   end_timing() { samples += MPI_Wtime() - t1; count ++; };

	double get_sum() { return samples; };
	double get_mean() { return get_sum() / count; };
	
	char* to_string() {
		sprintf(out_str, "%s: t(%lf) m(%lf)", name, get_sum(), get_mean());
		return out_str;
	};
	double get_low_mean_percentile(int percentile);
};

#endif /* FTLIB_BENCHMARK_HPP */
