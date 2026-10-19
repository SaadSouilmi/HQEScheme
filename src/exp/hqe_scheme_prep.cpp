#include <cstdio>
#include <complex.h>


#include "../base.h"

#include "../profile.h"
#include "../mittag_leffler.h"
#include "../hqe.h"

#include "../profile.cpp"
#include "../mittag_leffler.cpp"
#include "../hqe.cpp"



inline void KeepAlive(F64 X)
{
	__asm__ volatile("": : "x"(X));
}


I32
main()
{
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
	U64 NSteps = 1024;
	F64 Delta = 1.0 / NSteps;

	constexpr U64 Warmup = 1000;
	for (U64 J = 0; J < Warmup; J++)
	{
		KeepAlive(KernelEval(J*Delta, &Model));
	}

	constexpr U64 NIters   = 10000;

	profile_stats Stats = {};
	F64 CheckSum = 0;
	for (U64 J = 0; J < NIters; J++)
	{
		U64 Start = ReadTSC();
		hqe_scheme HQE;
		I32 Code = HQEInit(&HQE, 1.0, 1024, &Model);
		if (Code == -1) 
		{ 
			fprintf(stderr, "Malloc error\n");
			return -1;
		}
		U64 Elapsed = ReadTSC() - Start;
		IncrementStats(&Stats, Elapsed);
		CheckSum += HQE.CEps;
	}


	printf("CheckSum = %f\n\n", CheckSum);

	printf("Rough Hestom Model Setup: %.3f ± %.3f us/model\n", 
		TSCToWall(Stats.AvgTicks, TimeUnit_Us),
		TSCToWall(Stats.StdTicks, TimeUnit_Us));
}
