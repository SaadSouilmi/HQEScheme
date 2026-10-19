#include <cstdio>

#include "../base.h"
#include "../rng.h"
#include "../mittag_leffler.h"
#include "../hqe.h"
#include "../profile.h"

#include "../rng.cpp"
#include "../mittag_leffler.cpp"
#include "../hqe.cpp"
#include "../profile.cpp"

inline void KeepAlive(F64 X)
{
	__asm__ volatile("": : "x"(X));
}

I32
main()
{
	constexpr U64 Warmup   = 100;
	constexpr U64 NPaths   = 10000;
	
	model_params Params = {
		.Alpha=0.55,
		.V0=0.4,
		.Theta=0.2,
		.Lambda=0.33,
		.Nu=0.1,
		.Rho=-0.8
	};
	rough_heston_model Model;
	ModelInit(&Model, &Params);
	constexpr U64 NSteps = 1024;
	F64 Delta = 1.0 / NSteps;
	hqe_scheme HQE;
	HQEInit(&HQE, 1.0, NSteps, &Model);
	rng Rng;
	RandomSeed(&Rng, 1337, 0);
	
	F64 CheckSum = 0.0;
	for (U64 J = 0; J < Warmup; J++)
	{
		F64 *Var  = (F64*)malloc(sizeof(F64)*(NSteps+1));
		F64 *LogS = (F64*)malloc(sizeof(F64)*(NSteps+1));
		I32 Status = SamplePath(Var, LogS, 100.0, &HQE, &Rng);
		if (Status == -1) { return -1; }
		free(Var); free(LogS);
		CheckSum += Var[NSteps];
	}

	profile_stats Stats = {};
	for (U64 J = 0; J < NPaths; J++)
	{
		U64 Start = ReadTSC();
		F64 *Var  = (F64*)malloc(sizeof(F64)*(NSteps+1));
		F64 *LogS = (F64*)malloc(sizeof(F64)*(NSteps+1));
		I32 Status = SamplePath(Var, LogS, 100.0, &HQE, &Rng);
		if (Status == -1) { return -1; }
		free(Var); free(LogS);
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&Stats, Elapsed);
		CheckSum += Var[NSteps];
	}

	KeepAlive(CheckSum);

	printf("HQE scheme: %.3f ± %.3f us/path, 10k paths : %.3f ms\n", 
		TSCToWall(Stats.AvgTicks, TimeUnit_Us),
		TSCToWall(Stats.StdTicks, TimeUnit_Us),
		TSCToWall(Stats.AvgTicks*NPaths, TimeUnit_Ms));


}
