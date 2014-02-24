#ifndef __TIMER_H__
#define __TIMER_H__

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#elif WIN64
#include <windows.h>
#include <stdio.h>
#else
#include <sys/time.h>
#endif

#ifdef WIN32
double PCFreq = 0.0;
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

inline
	void StartTimer()
{
#ifdef WIN32
	LARGE_INTEGER li;
	if(!QueryPerformanceFrequency(&li))
	printf("QueryPerformanceFrequency failed!\n");

	PCFreq = (double)li.QuadPart/1000.0;

	QueryPerformanceCounter(&li);
	timerStart = li.QuadPart;
#else
	gettimeofday(&timerStart, NULL);
#endif
}

// time elapsed in ms
inline
	double GetTimer()
{
#ifdef WIN32
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return (double)(li.QuadPart-timerStart)/PCFreq;
#else
	struct timeval timerStop, timerElapsed;
	gettimeofday(&timerStop, NULL);
	timersub(&timerStop, &timerStart, &timerElapsed);
	return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
#endif
}

#endif // __TIMER_H__