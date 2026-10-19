#ifndef HQE_H
#define HQE_H

struct model_params {
	F64 Alpha; // Alpha = H + 0.5
	F64 V0;
	F64 Theta;
	F64 Lambda;
	F64 Nu;
	F64 Rho;
};

struct rough_heston_model {
	model_params Params;
	MLQuad Q;
	MLQuad Q2;
};

struct kernel_sqr_int {
	const rough_heston_model *Model;
	F64 P; F64 L; F64 D; F64 C;
};

struct quad_result {
	F64 A; F64 B;
	F64 Value;
	F64 Error;
};

struct hqe_scheme {
	const rough_heston_model *Model;
	F64 Delta; U64 NSteps;
	F64 K1; F64 *W; F64 *Xi0;
	F64 Beta; F64 CChi; F64 CEps;
};


internal void ModelInit(rough_heston_model *Model, const model_params *Params);
internal inline F64 KernelEval(F64 t, const rough_heston_model *Model);
internal inline F64 KernelInt(F64 Delta, const rough_heston_model *Model);
internal void KernelSqrIntInit(kernel_sqr_int *KerSqr, F64 Delta, U64 J, const rough_heston_model *Model);
internal inline F64 KernelSqrIntEval(F64 U, const kernel_sqr_int *KerSqr);
internal quad_result Quadrature(F64 A, F64 B, const kernel_sqr_int *KerSqr);
internal F64 KernelSqrInt(const kernel_sqr_int *KerSqr);
internal inline F64 FwdVar0(F64 t, const rough_heston_model *Model);
internal I32 HQEInit(hqe_scheme *HQE, F64 T, U64 NSteps, const rough_heston_model *Model);

#endif // HQE_H
