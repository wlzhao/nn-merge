#ifndef TIMER_H
#define TIMER_H

#include <ctime>

class Timer
{
public:
	Timer();
    void start();
	tm  *end();
	void end(bool TXT_OUT);

    static void test();
	virtual ~Timer();
private:
	bool started;
	time_t start_t, end_t;
};

#endif
