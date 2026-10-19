#include <math.h>
#include <stdlib.h>

internal void
ModelInit(rough_heston_model *Model, const model_params *Params)
{
	Model->Params = *Params;
	MLQuadInit(Params->Alpha, Params->Alpha, &Model->Q);
	MLQuadInit(Params->Alpha, Params->Alpha+1, &Model->Q2);
}

internal inline F64
KernelEval(F64 t, const rough_heston_model *Model)
{
	if (Model->Params.Nu == 0) { return 0; }

	F64 ML = MLQuadEval(
		Model->Params.Lambda * pow(t, Model->Params.Alpha),
		&Model->Q
	);
	return Model->Params.Nu * pow(t, Model->Params.Alpha - 1) * ML;
}

// Comptue K1 = \int_0^\Delta kernel(t)dt
internal inline F64
KernelInt(F64 Delta, const rough_heston_model *Model)
{
	F64 DeltaPow = pow(Delta, Model->Params.Alpha);
	if (Model->Params.Lambda == 0)
	{
		return Model->Params.Nu 
		* DeltaPow
		* Model->Q2.AtZero;
	}

	F64 ML = MLQuadEval(
		Model->Params.Lambda * DeltaPow,
		&Model->Q2
	);

	return Model->Params.Nu * DeltaPow * ML;
}

internal void
KernelSqrIntInit(kernel_sqr_int *KerSqr, F64 Delta, U64 J, const rough_heston_model *Model)
{
	KerSqr->Model = Model;
	F64 JDelta = (F64)J * Delta;
	KerSqr->P = 2*Model->Params.Alpha - 1;
	KerSqr->L = pow(JDelta, KerSqr->P);
	KerSqr->D = J == 0
    ? pow(Delta, KerSqr->P)
		: KerSqr->L * expm1(KerSqr->P * log1p(1.0/(F64)J));
	KerSqr->C = Model->Params.Nu * Model->Params.Nu * KerSqr->D / KerSqr->P;
}

internal inline F64
KernelSqrIntEval(F64 U, const kernel_sqr_int *KerSqr)
{
	F64 ML = MLQuadEval(
		KerSqr->Model->Params.Lambda 
		* pow(KerSqr->L + KerSqr->D*U, KerSqr->Model->Params.Alpha/KerSqr->P),
		&KerSqr->Model->Q
	);

	return KerSqr->C * ML * ML;
}


internal quad_result
Quadrature(F64 A, F64 B, const kernel_sqr_int *KerSqr)
{
	read_only F64 Nodes[] = {0.9914553711208126, 0.9491079123427585, 0.8648644233597691,
		0.7415311855993945, 0.5860872354676911, 0.4058451513773972, 0.2077849550078985};
	read_only F64 KW[] = {0.02293532201052922, 0.06309209262997855, 0.1047900103222502,
		0.1406532597155259, 0.1690047266392679, 0.1903505780647854, 0.2044329400752989, 0.2094821410847278};
	read_only F64 GW[] = {0.1294849661688697, 0.2797053914892767, 0.3818300505051189, 0.4179591836734694};
	F64 Mid = (A + B)/2, Half = (B - A)/2, Center = KernelSqrIntEval(Mid, KerSqr);
	F64 K = KW[7]*Center, G = GW[3]*Center;
	for (I32 J = 0; J < 7; ++J)
	{
		F64 Pair = KernelSqrIntEval(Mid - Half*Nodes[J], KerSqr) 
				   + KernelSqrIntEval(Mid + Half*Nodes[J], KerSqr);
		K += KW[J]*Pair;
		if (J % 2) { G += GW[J/2]*Pair; }
	}
	K *= Half; G *= Half;
	return (quad_result){A, B, K, fabs(K - G)};
}

internal F64
KernelSqrInt(const kernel_sqr_int *KerSqr)
{
	constexpr U64 MaxSplits = 512;
	constexpr F64 AbsTol = 2e-13;
	constexpr F64 RelTol = 2e-12;

	quad_result Results[MaxSplits];
	Results[0] = Quadrature(0.0, 1.0, KerSqr);
	F64 Value = Results[0].Value;
	F64 Error = Results[0].Error;
	for (U64 Count = 1; ; Count++)
	{
		if (!isfinite(Value) || !isfinite(Error) || Count == MaxSplits) { return NAN; }
		if (Error <= AbsTol + 2e-12*fabs(Value)) { return Value; }
		U64 Worst = 0;
		for (U64 J = 1; J < Count; J++)
		{
			if (Results[J].Error > Results[Worst].Error)
			{
				Worst = J;
			}
		}
		quad_result Curr = Results[Worst];
		F64 Mid = (Curr.A + Curr.B)/2;
		if (Mid == Curr.A || Mid == Curr.B) { return NAN; }
		Results[Worst] = Quadrature(Curr.A, Mid, KerSqr);
		Results[Count] = Quadrature(Mid, Curr.B, KerSqr);
		Value += Results[Worst].Value + Results[Count].Value - Curr.Value;
		Error = Max(0, Error + Results[Worst].Error + Results[Count].Error - Curr.Error);
	}
	return NAN;
}

internal inline F64
FwdVar0(F64 t, const rough_heston_model *Model)
{
	F64 tpow = pow(t, Model->Params.Alpha);
	F64 ML = MLQuadEval(
		Model->Params.Lambda * tpow,
		&Model->Q2
	);

	return Model->Params.V0 + Model->Params.Lambda * (Model->Params.Theta - Model->Params.V0) * tpow * ML;
}

internal I32
HQEInit(hqe_scheme *HQE, F64 T, U64 NSteps, const rough_heston_model *Model)
{
	HQE->Model = Model;
	HQE->W = (F64*)malloc(sizeof(F64)*NSteps);
	if (!HQE->W) { return -1; }
	HQE->Xi0 = (F64*)malloc(sizeof(F64)*(NSteps+1));
	if (!HQE->Xi0) { HQE->W = nullptr; free(HQE->W); return -1; }

	F64 Delta = T / (F64)NSteps;
	HQE->Delta = Delta;
	HQE->NSteps = NSteps;
	
	HQE->K1 = KernelInt(Delta, Model);
	for (U64 J = 0; J < NSteps; J++)
	{
		kernel_sqr_int KerSqr = {};
		KernelSqrIntInit(&KerSqr, Delta, J, Model);
		F64 K2 = KernelSqrInt(&KerSqr);
		if (J == 0) { HQE->CEps = K2; }
		HQE->W[J] = sqrt(K2 / Delta);
		HQE->Xi0[J] = FwdVar0(J*Delta, Model);
	}
	HQE->Xi0[NSteps] = FwdVar0(T, Model);
	HQE->Beta = HQE->K1 / Delta;
	HQE->CChi = HQE->K1 * HQE->K1 / Delta;
	HQE->CEps -= HQE->CChi;
	return 0;
}

