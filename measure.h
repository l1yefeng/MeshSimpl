// Created by nickl (LYF) on 12/21/18.
//
// measure.h defines a class Measure helping measure execution time
//
// Example:
// 	measure.start();
// 	do_task1();
// 	long ms1 = measure.stop();
// 	do_task2();
// 	long ms2 = measure.stop();
// 	cout << ms1 << "ms doing task1; " << ms2 << "ms doing task2"

#ifndef MESH_SIMPL_MEASURE_H
#define MESH_SIMPL_MEASURE_H

#include <chrono>

class Measure
{
public:
	Measure() :before(std::chrono::steady_clock::now()) { }
	void start() { before = std::chrono::steady_clock::now(); }
	long stop()
	{
		const auto after = std::chrono::steady_clock::now();
		const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(after-before);
		before = after;
		return ms.count();
	}
private:
	std::chrono::steady_clock::time_point before;
};

#endif // MESH_SIMPL_MEASURE_H
