#ifndef PROFILE_H
#define PROFILE_H

enum time_unit {
	TimeUnit_Sec,
	TimeUnit_Ms,
	TimeUnit_Us,
	TimeUnit_Ns,
};

internal const char *
TimeUnitName(time_unit Unit)
{
    switch (Unit)
    {
        case TimeUnit_Sec:  return "s";
        case TimeUnit_Ms:   return "ms";
        case TimeUnit_Us:   return "us";
        case TimeUnit_Ns:   return "ns";
    }

    return "?";
}

struct tsc_calibration {
	U64 TSCTicks;
	U64 OSNanos;
};

struct profile_stats {
    U64 Count;

    F64 AvgTicks;
    F64 M2Ticks;
    F64 StdTicks;

    U64 MinTicks;
    U64 MaxTicks;
};

internal U64 ReadTSC();
internal U64 ReadOSTime();
internal tsc_calibration GetTSCCalibration();
internal F64 TSCToWall(F64 ElapsedTicks, time_unit Unit);
internal F64 TSCToWall(U64 ElapsedTicks, time_unit Unit);
internal void IncrementStats(profile_stats *Stats, U64 Elapsed);

#endif // PROFILE_H


