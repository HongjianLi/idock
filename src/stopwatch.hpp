#pragma once
#ifndef IDOCK_STOPWATCH_HPP
#define IDOCK_STOPWATCH_HPP

#include <chrono>
#include <string>
using namespace std;
using namespace std::chrono;

//! Represents a stopwatch.
class stopwatch
{
public:
	//! Create a new stopwatch and start running immediately.
	static stopwatch start_new();

public:
	//! Constructs a standby stopwatch.
	explicit stopwatch();

	void start();

	void stop();

	void restart();

	void reset();

	bool is_running() const;

	//! Accumulated elapsed nanoseconds.
	long long elapsed() const;

	//! Accumulated elapsed seconds.
	double elapsed_sec() const;

	//! Accumulated elapsed time in mm:ss.sss format.
	string elapsed_str() const;

private:
	using clock = high_resolution_clock;
	size_t elapsed_ns; //!< Accumulated elapsed nanoseconds.
	time_point<clock> started_time; //!< Started time of the current running.
	bool running; //!< Current running status.
};

#endif
