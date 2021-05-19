#include "stopwatch.hpp"

//! Constructs a standby stopwatch.
stopwatch::stopwatch()
	: elapsed_ns(0)
	, started_time()
	, running(false)
{
}

stopwatch stopwatch::start_new()
{
	stopwatch sw;
	sw.start();
	return sw;
}

void stopwatch::start()
{
	if (running)
		return;
	running = true;
	started_time = clock::now();
}

void stopwatch::stop()
{
	if (!running)
		return;
	running = false;
	elapsed_ns += (clock::now() - started_time).count();
}

void stopwatch::restart()
{
	running = true;
	elapsed_ns = 0;
	started_time = clock::now();
}

void stopwatch::reset()
{
	running = false;
	elapsed_ns = 0;
}

bool stopwatch::is_running() const
{
	return running;
}

//! Accumulated elapsed nanoseconds.
long long stopwatch::elapsed() const
{
	if (running)
		return elapsed_ns + (clock::now() - started_time).count();
	return elapsed_ns;
}

//! Accumulated elapsed seconds.
double stopwatch::elapsed_sec() const
{
	return elapsed() / 1e9;
}

//! Accumulated elapsed time in mm:ss.sss format.
string stopwatch::elapsed_str() const
{
	auto e = elapsed(), sec = e / nanoseconds::period::den;
	auto fraction = e % nanoseconds::period::den / micro::den;
	auto min = sec / minutes::period::num;
	return (min < 10 ? "0" : "") + to_string(min) + ':' + (sec < 10 ? "0" : "") + to_string(sec) + to_string(fraction / 1000.0).substr(1);
}