#include <stdint.h>
#include <sys/time.h>
#include"TimeUtility.h"

struct timeval  start;

void TimeUtility::StartCounterMicro()
{
	gettimeofday(&start, NULL);
}
double TimeUtility::GetCounterMicro()
{
	struct timeval  end;
	gettimeofday(&end, NULL);

	long seconds  = end.tv_sec  - start.tv_sec;
	long useconds = end.tv_usec - start.tv_usec;
	long mtime = ((seconds) * 1000000 + useconds) + 0.5;

	return (double)mtime;
}

void TimeUtility::StartCounterMill()
{
	gettimeofday(&start, NULL);
}

double TimeUtility::GetCounterMill()
{
	struct timeval  end;
	gettimeofday(&end, NULL);

	long seconds  = end.tv_sec  - start.tv_sec;
	long useconds = end.tv_usec - start.tv_usec;
	long mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    return (double)mtime;
}
