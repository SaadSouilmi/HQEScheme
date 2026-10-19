// ======================================================================================
// Copyright © 2013-2014 Peter Jäckel.
//
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
// ======================================================================================

#ifndef JACKEL_H
#define JACKEL_H

#include <float.h>
#include <math.h>

enum option_type
{
	Option_Call = 1,
	Option_Put = -1,
};

#define JackelOneOverSqrtTwo 0.7071067811865475244008443621048490392848359376887
#define JackelOneOverSqrtTwoPi 0.3989422804014326779399460599343818684758586311649
#define JackelSqrtTwoPi 2.506628274631000502415765284811045253006986740610
#define JackelTwoPi 6.283185307179586476925286766559005768394338798750
#define JackelSqrtPiOverTwo 1.253314137315500251207882642405522626503493370305
#define JackelSqrtThree 1.732050807568877293527446341505872366942805253810
#define JackelSqrtOneOverThree 0.577350269189625764509148780501957455647601751270
#define JackelTwoPiOverSqrtTwentySeven 1.209199576156145233729385505094770488189377498728
#define JackelPiOverSix 0.523598775598298873077107230546583814032861566563

#define JackelSqrtF64Epsilon 1.4901161193847656e-08
#define JackelFourthRootF64Epsilon 0.0001220703125
#define JackelEighthRootF64Epsilon 0.011048543456039806
#define JackelSixteenthRootF64Epsilon 0.10511205190671433
#define JackelSqrtF64Min 1.4916681462400413e-154
#define JackelSqrtF64Max 1.3407807929942596e+154

#define JackelNormCdfAsymptoticExpansionFirstThreshold -10.0
#define JackelNormCdfAsymptoticExpansionSecondThreshold -67108864.0
#define JackelMinimumRationalCubicControlParameterValue -0.99999998509883881
#define JackelMaximumRationalCubicControlParameterValue 4.0564819207303341e+31
#define JackelAsymptoticExpansionAccuracyThreshold -10.0
#define JackelSmallTExpansionThreshold (2.0 * JackelSixteenthRootF64Epsilon)

internal F64
JackelCodyErfc(F64 Z);

internal F64
JackelCodyErfcx(F64 Z);

internal F64
JackelNormalCdf(F64 Z);

internal F64
JackelInverseNormalCdf(F64 U);

internal inline F64
JackelSquare(F64 X)
{
	return X * X;
}

internal inline F64
JackelNormalPdf(F64 X)
{
	return JackelOneOverSqrtTwoPi * exp(-.5 * X * X);
}

internal inline I32
JackelIsZero(F64 X)
{
	return fabs(X) < DBL_MIN;
}

internal inline F64
JackelHouseholderFactor(F64 Newton, F64 Halley, F64 Hh3)
{
	return (1 + 0.5 * Halley * Newton) / (1 + Newton * (Halley + Hh3 * Newton / 6));
}

internal F64
JackelRationalCubic(F64 X, F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 R);

internal F64
JackelMinimumCubicControl(F64 DL, F64 DR, F64 S, I32 PreferShapePreservation);

internal F64
JackelConvexCubicControlLeft(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeL,
	I32 PreferShapePreservation);

internal F64
JackelConvexCubicControlRight(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeR,
	I32 PreferShapePreservation);

internal F64
JackelNormalisedIntrinsic(F64 X, F64 Q);

internal F64
JackelNormalisedBlackCall(F64 X, F64 S);

internal F64
JackelNormalisedVega(F64 X, F64 S);

internal void
JackelLowerMap(F64 X, F64 S, F64 *F, F64 *Fp, F64 *Fpp);

internal F64
JackelInverseLowerMap(F64 X, F64 F);

internal void
JackelUpperMap(F64 X, F64 S, F64 *F, F64 *Fp, F64 *Fpp);

internal F64
JackelInverseUpperMap(F64 F);

internal F64
JackelImpliedVolatility(F64 Price, F64 Forward, F64 Strike, F64 Expiry, option_type Type);

internal void
JackelImpliedVolVec(const F64 *Prices, const F64 *Strikes, const option_type *Types,
	F64 Forward, F64 Expiry, F64 *Out, U64 N);

#endif // JACKEL_H

