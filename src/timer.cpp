#include "timer.h"

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

Timer::Timer()
{
    this->started = false;
}

void Timer::start()
{
    this->started = true;
    time(&this->start_t);
    localtime(&this->start_t);
}

tm* Timer::end()
{
    tm *tm_s, *tm_e, *tm_elaps=(tm *)malloc(sizeof(tm));
    int min = 0, sec = 0, hour = 0, yday = 0;
    tm_elaps->tm_hour  = 0;
    tm_elaps->tm_isdst = 0;
    tm_elaps->tm_mday  = 0;
    tm_elaps->tm_min   = 0;
    tm_elaps->tm_mon   = 0;
    tm_elaps->tm_sec   = 0;
    tm_elaps->tm_wday  = 0;
    tm_elaps->tm_yday  = 0;
    tm_elaps->tm_year  = 0;

    if(this->started)
    {
        time(&this->end_t);
        tm_s = localtime(&this->start_t);
        yday = tm_s->tm_yday;
        hour = tm_s->tm_hour;
        min  = tm_s->tm_min;
        sec  = tm_s->tm_sec;
        tm_e = localtime(&this->end_t);
        tm_elaps->tm_yday = tm_e->tm_yday - yday;
        tm_elaps->tm_hour = tm_e->tm_hour - hour;
        tm_elaps->tm_min  = tm_e->tm_min  - min;
        tm_elaps->tm_sec  = tm_e->tm_sec  - sec;

        if(tm_elaps->tm_sec < 0)
        {
            tm_elaps->tm_min = tm_elaps->tm_min - 1;
            tm_elaps->tm_sec = tm_elaps->tm_sec+60;
        }

        if(tm_elaps->tm_min < 0)
        {
            tm_elaps->tm_hour = tm_elaps->tm_hour-1;
            tm_elaps->tm_min  = tm_elaps->tm_min+60;
        }

        if(tm_elaps->tm_hour < 0)
        {
            tm_elaps->tm_yday = tm_elaps->tm_yday - 1;
            tm_elaps->tm_hour = tm_elaps->tm_hour + 24;
        }
    }
    return tm_elaps;
}

void Timer::end(bool TXT_OUT)
{
    tm *tm_s,*tm_e,*tm_elaps = (tm *)malloc(sizeof(tm));
    int min = 0, sec = 0, hour = 0, yday = 0;
    tm_elaps->tm_hour  = 0;
    tm_elaps->tm_isdst = 0;
    tm_elaps->tm_mday  = 0;
    tm_elaps->tm_min   = 0;
    tm_elaps->tm_mon   = 0;
    tm_elaps->tm_sec   = 0;
    tm_elaps->tm_wday  = 0;
    tm_elaps->tm_yday  = 0;
    tm_elaps->tm_year  = 0;

    if(this->started)
    {
        time(&this->end_t);
        tm_s = localtime(&this->start_t);
        yday = tm_s->tm_yday;
        hour = tm_s->tm_hour;
        min  = tm_s->tm_min;
        sec  = tm_s->tm_sec;
        tm_e = localtime(&this->end_t);
        tm_elaps->tm_yday = tm_e->tm_yday-yday;
        tm_elaps->tm_hour = tm_e->tm_hour-hour;
        tm_elaps->tm_min  = tm_e->tm_min-min;
        tm_elaps->tm_sec  = tm_e->tm_sec-sec;

        if(tm_elaps->tm_sec < 0)
        {
            tm_elaps->tm_min = tm_elaps->tm_min-1;
            tm_elaps->tm_sec = tm_elaps->tm_sec+60;
        }

        if(tm_elaps->tm_min < 0)
        {
            tm_elaps->tm_hour = tm_elaps->tm_hour-1;
            tm_elaps->tm_min  = tm_elaps->tm_min+60;
        }

        if(tm_elaps->tm_hour < 0)
        {
            tm_elaps->tm_yday = tm_elaps->tm_yday - 1;
            tm_elaps->tm_hour = tm_elaps->tm_hour + 24;
        }
    }
    if(TXT_OUT)
    {
        cout<<endl<<"Time costs:";
        if(tm_elaps->tm_yday != 0.0f)
        {
            cout<<" "<<tm_elaps->tm_yday<<" day(s)";
        }

        if(tm_elaps->tm_hour != 0.0f)
        {
            cout<<"  "<<tm_elaps->tm_hour<<" hour(s)";
        }
        if(tm_elaps->tm_min != 0.0f)
        {
            cout<<"  "<< tm_elaps->tm_min<<" minute(s)";
        }
        if(tm_elaps->tm_sec != 0.0f)
        {
            cout<<"  "<<tm_elaps->tm_sec<<" second(s)";
        }
        cout<<endl;
    }
    delete tm_elaps;
    return ;
}

void Timer::test()
{
    Timer *my = new Timer();
    my->start();
    for(int i = 0; i < 100; i++)
        for(int j = 0; j < 3000; j++)
        {
            cout<<"Hellow"<<endl;
        }
    my->end(true);
    return ;
}

Timer::~Timer()
{

}
