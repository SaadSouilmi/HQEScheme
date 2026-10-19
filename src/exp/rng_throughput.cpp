#include "../base.h"

#include "../rng.h"
#include "../profile.h"

#include "../rng.cpp"
#include "../profile.cpp"

#include <cstdio>


inline void KeepAlive(F64 X)
{
	__asm__ volatile("": : "x"(X));
}

inline void KeepAlive(U64 X)
{
	__asm__ volatile("": : "x"(X));
}


I32
main()
{
	rng Rng;
	RandomSeed(&Rng, 1337, 0);

	GetTSCCalibration();
	constexpr U64 Warmup = 100000;
	for (U64 i = 0; i < Warmup; i++)
	{
		KeepAlive(RandomNormal(&Rng));
		KeepAlive(RandomExponential(&Rng));
	}

	constexpr U64 NBatches = 10000;
	constexpr U64 NSamples = 10000;

	profile_stats U64Profile = {};
	for (U64 B = 0; B < NBatches; B++)
	{
		U64 Start = ReadTSC();
		for (U64 t = 0; t < NSamples; t++)
		{
			KeepAlive(RandomU64(&Rng));
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&U64Profile, Elapsed);
	}
	
	printf("U64 throughput: %.3f samples/sec\n", (F64)NSamples/TSCToWall(U64Profile.AvgTicks, TimeUnit_Sec));
	printf("U64 avg sampling cost: %.3f ± %.3f ns/sample\n\n", TSCToWall(U64Profile.AvgTicks, TimeUnit_Ns)/(F64)NSamples, TSCToWall(U64Profile.StdTicks, TimeUnit_Ns)/(F64)NSamples);

	profile_stats UniformProfile = {};
	for (U64 B = 0; B < NBatches; B++)
	{
		U64 Start = ReadTSC();
		for (U64 t = 0; t < NSamples; t++)
		{
			KeepAlive(RandomUniform(&Rng));
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&UniformProfile, Elapsed);
	}
	
	printf("Uniform throughput: %.3f samples/sec\n", (F64)NSamples/TSCToWall(UniformProfile.AvgTicks, TimeUnit_Sec));
	printf("Uniform avg sampling cost: %.3f ± %.3f ns/sample\n\n", TSCToWall(UniformProfile.AvgTicks, TimeUnit_Ns)/(F64)NSamples, TSCToWall(UniformProfile.StdTicks, TimeUnit_Ns)/(F64)NSamples);

	profile_stats NormProfile = {};
	for (U64 B = 0; B < NBatches; B++)
	{
		U64 Start = ReadTSC();
		for (U64 t = 0; t < NSamples; t++)
		{
			KeepAlive(RandomNormal(&Rng));
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&NormProfile, Elapsed);
	}
	
	printf("Normal throughput: %.3f samples/sec\n", (F64)NSamples/TSCToWall(NormProfile.AvgTicks, TimeUnit_Sec));
	printf("Normal avg sampling cost: %.3f ± %.3f ns/sample\n\n", TSCToWall(NormProfile.AvgTicks, TimeUnit_Ns)/(F64)NSamples, TSCToWall(NormProfile.StdTicks, TimeUnit_Ns)/(F64)NSamples);

	profile_stats ExpProfile = {};
	for (U64 B = 0; B < NBatches; B++)
	{
		U64 Start = ReadTSC();
		for (U64 t = 0; t < NSamples; t++)
		{
			KeepAlive(RandomExponential(&Rng));
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&ExpProfile, Elapsed);
	}
	
	printf("Exp throughput: %.3f samples/sec\n", (F64)NSamples/TSCToWall(ExpProfile.AvgTicks, TimeUnit_Sec));
	printf("Exp avg sampling cost: %.3f ± %.3f ns/sample\n\n", TSCToWall(ExpProfile.AvgTicks, TimeUnit_Ns)/(F64)NSamples, TSCToWall(ExpProfile.StdTicks, TimeUnit_Ns)/(F64)NSamples);

	std::pair<F64, F64> Psi[NSamples];
	for (U64 J = 0; J < NSamples; J++)
	{
		bool Flag = RandomUniform(&Rng) < 0.5;
		Psi[J].first = (Flag) ? 1 : 2;
		Psi[J].second = 2;
	}

	profile_stats QEProfile = {};
	for (U64 B = 0; B < NBatches; B++)
	{
		U64 Start = ReadTSC();
		for (U64 t = 0; t < NSamples; t++)
		{
			KeepAlive(RandomQE(Psi[t].first, Psi[t].second, &Rng));
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&QEProfile, Elapsed);
	}

	printf("QE throughput: %.3f samples/sec\n", (F64)NSamples/TSCToWall(QEProfile.AvgTicks, TimeUnit_Sec));
	printf("QE avg sampling cost: %.3f ± %.3f ns/sample\n\n", 
		   TSCToWall(QEProfile.AvgTicks, TimeUnit_Ns)/(F64)NSamples, 
		   TSCToWall(QEProfile.StdTicks, TimeUnit_Ns)/(F64)NSamples);
}
