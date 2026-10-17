// ======================================================================================
// Copyright © 2013-2014 Peter Jäckel.
//
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
//
// WARRANTY DISCLAIMER
// The Software is provided "as is" without warranty of any kind, either express
// or implied, including without limitation any implied warranties of condition,
// uninterrupted use, merchantability, fitness for a particular purpose, or
// non-infringement.
// ======================================================================================

#include "jackel.h"

// W. J. Cody, Argonne National Laboratory, March 19, 1990.
// Original Fortran: http://www.netlib.org/specfun/erf
// Translated by f2c (version 20100827) and adapted by hand.
// Rational Chebyshev approximations for the error function,
// Mathematics of Computation (1969), pp. 631-638.

internal F64
JackelCodyErf(F64 X, const I32 Mode)
{
	read_only F64 A[5] =
	{
		3.1611237438705656, 113.864154151050156, 377.485237685302021,
		3209.37758913846947, .185777706184603153
	};
	read_only F64 B[4] =
	{
		23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062
	};
	read_only F64 C[9] =
	{
		.564188496988670089, 8.88314979438837594, 66.1191906371416295,
		298.635138197400131, 881.95222124176909, 1712.04761263407058, 2051.07837782607147,
		1230.33935479799725, 2.15311535474403846e-8
	};
	read_only F64 D[8] =
	{
		15.7449261107098347, 117.693950891312499, 537.181101862009858,
		1621.38957456669019, 3290.79923573345963, 4362.61909014324716, 3439.36767414372164,
		1230.33935480374942
	};
	read_only F64 P[6] =
	{
		.305326634961232344, .360344899949804439, .125781726111229246,
		.0160837851487422766, 6.58749161529837803e-4, .0163153871373020978
	};
	read_only F64 Q[5] =
	{
		2.56852019228982242, 1.87295284992346047, .527905102951428412,
		.0605183413124413191, .00233520497626869185
	};

	read_only F64 Zero = 0.;
	read_only F64 Half = .5;
	read_only F64 One = 1.;
	read_only F64 Two = 2.;
	read_only F64 Four = 4.;
	read_only F64 InvSqrtPi = 0.56418958354775628695;
	read_only F64 Threshold = .46875;
	read_only F64 Sixteen = 16.;

	F64 Y, Delta, YSquared, Denominator, Numerator, Result;

	read_only F64 XInf = 1.79e308;
	read_only F64 XNeg = -26.628;
	read_only F64 XSmall = 1.11e-16;
	read_only F64 XBig = 26.543;
	read_only F64 XHuge = 6.71e7;
	read_only F64 XMax = 2.53e307;
	Y = fabs(X);
	if (Y <= Threshold)
	{
		YSquared = Zero;
		if (Y > XSmall)
		{
			YSquared = Y * Y;
		}
		Numerator = A[4] * YSquared;
		Denominator = YSquared;
		for (I32 I = 1; I <= 3; ++I)
		{
			Numerator = (Numerator + A[I - 1]) * YSquared;
			Denominator = (Denominator + B[I - 1]) * YSquared;
		}
		Result = X * (Numerator + A[3]) / (Denominator + B[3]);
		if (Mode != 0)
		{
			Result = One - Result;
		}
		if (Mode == 2)
		{
			Result = exp(YSquared) * Result;
		}
		goto Done;
	}
	else if (Y <= Four)
	{
		Numerator = C[8] * Y;
		Denominator = Y;
		for (I32 I = 1; I <= 7; ++I)
		{
			Numerator = (Numerator + C[I - 1]) * Y;
			Denominator = (Denominator + D[I - 1]) * Y;
		}
		Result = (Numerator + C[7]) / (Denominator + D[7]);
		if (Mode != 2)
		{
			F64 Temp = Y * Sixteen;
			YSquared = trunc(Temp) / Sixteen;
			Delta = (Y - YSquared) * (Y + YSquared);
			Temp = exp(-YSquared * YSquared) * exp(-Delta);
			Result = Temp * Result;
		}
	}
	else
	{
		Result = Zero;
		if (Y >= XBig)
		{
			if (Mode != 2 || Y >= XMax)
			{
				goto FixSign;
			}
			if (Y >= XHuge)
			{
				Result = InvSqrtPi / Y;
				goto FixSign;
			}
		}
		YSquared = One / (Y * Y);
		Numerator = P[5] * YSquared;
		Denominator = YSquared;
		for (I32 I = 1; I <= 4; ++I)
		{
			Numerator = (Numerator + P[I - 1]) * YSquared;
			Denominator = (Denominator + Q[I - 1]) * YSquared;
		}
		Result = YSquared * (Numerator + P[4]) / (Denominator + Q[4]);
		Result = (InvSqrtPi - Result) / Y;
		if (Mode != 2)
		{
			F64 Temp = Y * Sixteen;
			YSquared = trunc(Temp) / Sixteen;
			Delta = (Y - YSquared) * (Y + YSquared);
			Temp = exp(-YSquared * YSquared) * exp(-Delta);
			Result = Temp * Result;
		}
	}
FixSign:
	if (Mode == 0)
	{
		Result = (Half - Result) + Half;
		if (X < Zero)
		{
			Result = -(Result);
		}
	}
	else if (Mode == 1)
	{
		if (X < Zero)
		{
			Result = Two - Result;
		}
	}
	else
	{
		if (X < Zero)
		{
			if (X < XNeg)
			{
				Result = XInf;
			}
			else
			{
				F64 Temp = X * Sixteen;
				YSquared = trunc(Temp) / Sixteen;
				Delta = (X - YSquared) * (X + YSquared);
				Y = exp(YSquared * YSquared) * exp(Delta);
				Result = Y + Y - Result;
			}
		}
	}
Done:
	return Result;
}

internal F64
JackelCodyErfc(F64 X)
{
	return JackelCodyErf(X, 1);
}

internal F64
JackelCodyErfcx(F64 X)
{
	return JackelCodyErf(X, 2);
}

// Abramowitz and Stegun, Pocketbook of Mathematical Functions,
// equation (26.2.12), p. 408, ISBN 3-87144818-4.
internal F64
JackelNormalCdf(F64 Z)
{
	if (Z <= JackelNormCdfAsymptoticExpansionFirstThreshold)
	{
		F64 Sum = 1;
		if (Z >= JackelNormCdfAsymptoticExpansionSecondThreshold)
		{
			F64 ZSquared = Z * Z, I = 1, G = 1, X, Y, A = DBL_MAX, LastA;
			do
			{
				LastA = A;
				X = (4 * I - 3) / ZSquared;
				Y = X * ((4 * I - 1) / ZSquared);
				A = G * (X - Y);
				Sum -= A;
				G *= Y;
				++I;
				A = fabs(A);
			} while (LastA > A && A >= fabs(Sum * DBL_EPSILON));
		}
		return -JackelNormalPdf(Z) * Sum / Z;
	}
	return 0.5 * JackelCodyErfc(-Z * JackelOneOverSqrtTwo);
}

// Algorithm AS241, Applied Statistics (1988), vol. 37, no. 3.
// http://lib.stat.cmu.edu/apstat/241
internal F64
JackelInverseNormalCdf(F64 U)
{
	const F64 Split1 = 0.425;
	const F64 Split2 = 5.0;
	const F64 Const1 = 0.180625;
	const F64 Const2 = 1.6;

	const F64 A0 = 3.3871328727963666080E0;
	const F64 A1 = 1.3314166789178437745E+2;
	const F64 A2 = 1.9715909503065514427E+3;
	const F64 A3 = 1.3731693765509461125E+4;
	const F64 A4 = 4.5921953931549871457E+4;
	const F64 A5 = 6.7265770927008700853E+4;
	const F64 A6 = 3.3430575583588128105E+4;
	const F64 A7 = 2.5090809287301226727E+3;
	const F64 B1 = 4.2313330701600911252E+1;
	const F64 B2 = 6.8718700749205790830E+2;
	const F64 B3 = 5.3941960214247511077E+3;
	const F64 B4 = 2.1213794301586595867E+4;
	const F64 B5 = 3.9307895800092710610E+4;
	const F64 B6 = 2.8729085735721942674E+4;
	const F64 B7 = 5.2264952788528545610E+3;
	const F64 C0 = 1.42343711074968357734E0;
	const F64 C1 = 4.63033784615654529590E0;
	const F64 C2 = 5.76949722146069140550E0;
	const F64 C3 = 3.64784832476320460504E0;
	const F64 C4 = 1.27045825245236838258E0;
	const F64 C5 = 2.41780725177450611770E-1;
	const F64 C6 = 2.27238449892691845833E-2;
	const F64 C7 = 7.74545014278341407640E-4;
	const F64 D1 = 2.05319162663775882187E0;
	const F64 D2 = 1.67638483018380384940E0;
	const F64 D3 = 6.89767334985100004550E-1;
	const F64 D4 = 1.48103976427480074590E-1;
	const F64 D5 = 1.51986665636164571966E-2;
	const F64 D6 = 5.47593808499534494600E-4;
	const F64 D7 = 1.05075007164441684324E-9;
	const F64 E0 = 6.65790464350110377720E0;
	const F64 E1 = 5.46378491116411436990E0;
	const F64 E2 = 1.78482653991729133580E0;
	const F64 E3 = 2.96560571828504891230E-1;
	const F64 E4 = 2.65321895265761230930E-2;
	const F64 E5 = 1.24266094738807843860E-3;
	const F64 E6 = 2.71155556874348757815E-5;
	const F64 E7 = 2.01033439929228813265E-7;
	const F64 F1 = 5.99832206555887937690E-1;
	const F64 F2 = 1.36929880922735805310E-1;
	const F64 F3 = 1.48753612908506148525E-2;
	const F64 F4 = 7.86869131145613259100E-4;
	const F64 F5 = 1.84631831751005468180E-5;
	const F64 F6 = 1.42151175831644588870E-7;
	const F64 F7 = 2.04426310338993978564E-15;

	if (U <= 0)
	{
		return log(U);
	}
	if (U >= 1)
	{
		return log(1 - U);
	}

	const F64 Q = U - 0.5;
	if (fabs(Q) <= Split1)
	{
		const F64 R = Const1 - Q * Q;
		return Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) * R + A2) * R + A1) * R + A0) /
			(((((((B7 * R + B6) * R + B5) * R + B4) * R + B3) * R + B2) * R + B1) * R + 1.0);
	}
	else
	{
		F64 R = Q < 0.0 ? U : 1.0 - U;
		R = sqrt(-log(R));
		F64 Ret;
		if (R < Split2)
		{
			R = R - Const2;
			Ret = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3) * R + C2) * R + C1) * R + C0) /
				(((((((D7 * R + D6) * R + D5) * R + D4) * R + D3) * R + D2) * R + D1) * R + 1.0);
		}
		else
		{
			R = R - Split2;
			Ret = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3) * R + E2) * R + E1) * R + E0) /
				(((((((F7 * R + F6) * R + F5) * R + F4) * R + F3) * R + F2) * R + F1) * R + 1.0);
		}
		return Q < 0.0 ? -Ret : Ret;
	}
}

internal F64
JackelRationalCubic(F64 X, F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 R)
{
	const F64 H = (XR - XL);
	if (fabs(H) <= 0)
	{
		return 0.5 * (YL + YR);
	}
	const F64 T = (X - XL) / H;
	if (!(R >= JackelMaximumRationalCubicControlParameterValue))
	{
		const F64 T = (X - XL) / H, OneMinusT = 1 - T, T2 = T * T,
				  OneMinusT2 = OneMinusT * OneMinusT;
		return (YR * T2 * T + (R * YR - H * DR) * T2 * OneMinusT +
				   (R * YL + H * DL) * T * OneMinusT2 + YL * OneMinusT2 * OneMinusT) /
			(1 + (R - 3) * T * OneMinusT);
	}
	return YR * T + YL * (1 - T);
}

internal F64
JackelCubicControlLeft(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeL)
{
	const F64 H = (XR - XL), Numerator = 0.5 * H * SecondDerivativeL + (DR - DL);
	if (JackelIsZero(Numerator))
	{
		return 0;
	}
	const F64 Denominator = (YR - YL) / H - DL;
	if (JackelIsZero(Denominator))
	{
		return Numerator > 0 ? JackelMaximumRationalCubicControlParameterValue
							 : JackelMinimumRationalCubicControlParameterValue;
	}
	return Numerator / Denominator;
}

internal F64
JackelCubicControlRight(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeR)
{
	const F64 H = (XR - XL), Numerator = 0.5 * H * SecondDerivativeR + (DR - DL);
	if (JackelIsZero(Numerator))
	{
		return 0;
	}
	const F64 Denominator = DR - (YR - YL) / H;
	if (JackelIsZero(Denominator))
	{
		return Numerator > 0 ? JackelMaximumRationalCubicControlParameterValue
							 : JackelMinimumRationalCubicControlParameterValue;
	}
	return Numerator / Denominator;
}

internal F64
JackelMinimumCubicControl(F64 DL, F64 DR, F64 S, I32 PreferShapePreservation)
{
	const I32 Monotonic = DL * S >= 0 && DR * S >= 0, Convex = DL <= S && S <= DR,
			  Concave = DL >= S && S >= DR;
	if (!Monotonic && !Convex && !Concave)
	{
		return JackelMinimumRationalCubicControlParameterValue;
	}
	const F64 DRMinusDL = DR - DL, DRMinusS = DR - S, SMinusDL = S - DL;
	F64 R1 = -DBL_MAX, R2 = R1;
	if (Monotonic)
	{
		if (!JackelIsZero(S))
		{
			R1 = (DR + DL) / S;
		}
		else if (PreferShapePreservation)
		{
			R1 = JackelMaximumRationalCubicControlParameterValue;
		}
	}
	if (Convex || Concave)
	{
		if (!(JackelIsZero(SMinusDL) || JackelIsZero(DRMinusS)))
		{
			R2 = Max(fabs(DRMinusDL / DRMinusS), fabs(DRMinusDL / SMinusDL));
		}
		else if (PreferShapePreservation)
		{
			R2 = JackelMaximumRationalCubicControlParameterValue;
		}
	}
	else if (Monotonic && PreferShapePreservation)
	{
		R2 = JackelMaximumRationalCubicControlParameterValue;
	}
	return Max(JackelMinimumRationalCubicControlParameterValue, Max(R1, R2));
}

internal F64
JackelConvexCubicControlLeft(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeL,
	I32 PreferShapePreservation)
{
	const F64 R = JackelCubicControlLeft(XL, XR, YL, YR, DL, DR, SecondDerivativeL);
	const F64 RMin =
		JackelMinimumCubicControl(DL, DR, (YR - YL) / (XR - XL), PreferShapePreservation);
	return Max(R, RMin);
}

internal F64
JackelConvexCubicControlRight(F64 XL, F64 XR, F64 YL, F64 YR, F64 DL, F64 DR, F64 SecondDerivativeR,
	I32 PreferShapePreservation)
{
	const F64 R = JackelCubicControlRight(XL, XR, YL, YR, DL, DR, SecondDerivativeR);
	const F64 RMin =
		JackelMinimumCubicControl(DL, DR, (YR - YL) / (XR - XL), PreferShapePreservation);
	return Max(R, RMin);
}

internal F64
JackelNormalisedIntrinsic(F64 X, F64 Q)
{
	if (Q * X <= 0)
	{
		return 0;
	}
	const F64 X2 = X * X;
	if (X2 < 98 * JackelFourthRootF64Epsilon)
	{
		return fabs(Max((Q < 0 ? -1 : 1) * X *
			(1 + X2 * ((1.0 / 24.0) + X2 * ((1.0 / 1920.0) +
			X2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * X2)))), 0.0));
	}
	const F64 BMax = exp(0.5 * X), OneOverBMax = 1 / BMax;
	return fabs(Max((Q < 0 ? -1 : 1) * (BMax - OneOverBMax), 0.));
}

internal F64
JackelBlackCallAsymptotic(F64 H, F64 T)
{
	const F64 E = (T / H) * (T / H), R = ((H + T) * (H - T)), Q = (H / R) * (H / R);
	const F64 AsymptoticExpansionSum =
		(2.0 + Q * (-6.0E0 - 2.0 * E
		+ 3.0 * Q * (1.0E1 + E * (2.0E1 + 2.0 * E)
		+ 5.0 * Q * (-1.4E1 + E * (-7.0E1 + E * (-4.2E1 - 2.0 * E))
		+ 7.0 * Q * (1.8E1 + E * (1.68E2 + E * (2.52E2 + E * (7.2E1 + 2.0 * E)))
		+ 9.0 * Q * (-2.2E1 + E * (-3.3E2 + E * (-9.24E2 + E * (-6.6E2 + E * (-1.1E2 - 2.0 *
			E))))
		+ 1.1E1 * Q * (2.6E1 + E * (5.72E2 + E * (2.574E3 + E * (3.432E3 + E * (1.43E3 + E *
			(1.56E2 + 2.0 * E)))))
		+ 1.3E1 * Q * (-3.0E1 + E * (-9.1E2 + E * (-6.006E3 + E * (-1.287E4 + E * (-1.001E4 +
			E * (-2.73E3 + E * (-2.1E2 - 2.0 * E))))))
		+ 1.5E1 * Q * (3.4E1 + E * (1.36E3 + E * (1.2376E4 + E * (3.8896E4 + E * (4.862E4 + E
			* (2.4752E4 + E * (4.76E3 + E * (2.72E2 + 2.0 * E)))))))
		+ 1.7E1 * Q * (-3.8E1 + E * (-1.938E3 + E * (-2.3256E4 + E * (-1.00776E5 + E *
			(-1.84756E5 + E * (-1.51164E5 + E * (-5.4264E4 + E * (-7.752E3 + E * (-3.42E2 - 2.0 *
			E))))))))
		+ 1.9E1 * Q * (4.2E1 + E * (2.66E3 + E * (4.0698E4 + E * (2.3256E5 + E * (5.8786E5 + E
			* (7.05432E5 + E * (4.0698E5 + E * (1.08528E5 + E * (1.197E4 + E * (4.2E2 + 2.0 *
			E)))))))))
		+ 2.1E1 * Q * (-4.6E1 + E * (-3.542E3 + E * (-6.7298E4 + E * (-4.90314E5 + E *
			(-1.63438E6 + E * (-2.704156E6 + E * (-2.288132E6 + E * (-9.80628E5 + E * (-2.01894E5
			+ E * (-1.771E4 + E * (-5.06E2 - 2.0 * E))))))))))
		+ 2.3E1 * Q * (5.0E1 + E * (4.6E3 + E * (1.0626E5 + E * (9.614E5 + E * (4.08595E6 + E
			* (8.9148E6 + E * (1.04006E7 + E * (6.53752E6 + E * (2.16315E6 + E * (3.542E5 + E *
			(2.53E4 + E * (6.0E2 + 2.0 * E)))))))))))
		+ 2.5E1 * Q * (-5.4E1 + E * (-5.85E3 + E * (-1.6146E5 + E * (-1.77606E6 + E *
			(-9.37365E6 + E * (-2.607579E7 + E * (-4.01166E7 + E * (-3.476772E7 + E *
			(-1.687257E7 + E * (-4.44015E6 + E * (-5.9202E5 + E * (-3.51E4 + E * (-7.02E2 - 2.0 *
			E))))))))))))
		+ 2.7E1 * Q * (5.8E1 + E * (7.308E3 + E * (2.3751E5 + E * (3.12156E6 + E * (2.003001E7
			+ E * (6.919458E7 + E * (1.3572783E8 + E * (1.5511752E8 + E * (1.0379187E8 + E *
			(4.006002E7 + E * (8.58429E6 + E * (9.5004E5 + E * (4.7502E4 + E * (8.12E2 + 2.0 *
			E)))))))))))))
		+ 2.9E1 * Q * (-6.2E1 + E * (-8.99E3 + E * (-3.39822E5 + E * (-5.25915E6 + E *
			(-4.032015E7 + E * (-1.6934463E8 + E * (-4.1250615E8 + E * (-6.0108039E8 + E *
			(-5.3036505E8 + E * (-2.8224105E8 + E * (-8.870433E7 + E * (-1.577745E7 + E *
			(-1.472562E6 + E * (-6.293E4 + E * (-9.3E2 - 2.0 * E))))))))))))))
		+ 3.1E1 * Q * (6.6E1 + E * (1.0912E4 + E * (4.74672E5 + E * (8.544096E6 + E *
			(7.71342E7 + E * (3.8707344E8 + E * (1.14633288E9 + E * (2.07431664E9 + E *
			(2.33360622E9 + E * (1.6376184E9 + E * (7.0963464E8 + E * (1.8512208E8 + E *
			(2.7768312E7 + E * (2.215136E6 + E * (8.184E4 + E * (1.056E3 + 2.0 * E)))))))))))))))
			+ 3.3E1 * (-7.0E1 + E * (-1.309E4 + E * (-6.49264E5 + E * (-1.344904E7 + E *
			(-1.4121492E8 + E * (-8.344518E8 + E * (-2.9526756E9 + E * (-6.49588632E9 + E *
			(-9.0751353E9 + E * (-8.1198579E9 + E * (-4.6399188E9 + E * (-1.6689036E9 + E *
			(-3.67158792E8 + E * (-4.707164E7 + E * (-3.24632E6 + E * (-1.0472E5 + E * (-1.19E3 -
			2.0 * E))))))))))))))))) * Q)))))))))))))))));
	const F64 B =
		JackelOneOverSqrtTwoPi * exp((-0.5 * (H * H + T * T))) * (T / R) * AsymptoticExpansionSum;
	return fabs(Max(B, 0.));
}

internal F64
JackelBlackCallErfcx(F64 H, F64 T)
{
	const F64 B = 0.5 * exp(-0.5 * (H * H + T * T)) *
		(JackelCodyErfcx(-JackelOneOverSqrtTwo * (H + T)) -
			JackelCodyErfcx(-JackelOneOverSqrtTwo * (H - T)));
	return fabs(Max(B, 0.0));
}

internal F64
JackelBlackCallSmallT(F64 H, F64 T)
{
	const F64 A = 1 + H * (0.5 * JackelSqrtTwoPi) * JackelCodyErfcx(-JackelOneOverSqrtTwo * H),
			  W = T * T, H2 = H * H;
	const F64 Expansion =
		2 * T * (A
		+ W * ((-1 + 3 * A + A * H2) / 6
		+ W * ((-7 + 15 * A + H2 * (-1 + 10 * A + A * H2)) / 120
		+ W * ((-57 + 105 * A + H2 * (-18 + 105 * A + H2 * (-1 + 21 * A + A * H2))) / 5040
		+ W * ((-561 + 945 * A + H2 * (-285 + 1260 * A + H2 * (-33 + 378 * A + H2 * (-1 + 36 *
			A + A * H2)))) / 362880
		+ W * ((-6555 + 10395 * A + H2 * (-4680 + 17325 * A + H2 * (-840 + 6930 * A + H2 *
			(-52 + 990 * A + H2 * (-1 + 55 * A + A * H2))))) / 39916800 + ((-89055 + 135135 * A +
			H2 * (-82845 + 270270 * A + H2 * (-20370 + 135135 * A + H2 * (-1926 + 25740 * A + H2
			* (-75 + 2145 * A + H2 * (-1 + 78 * A + A * H2)))))) * W) / 6227020800.0))))));
	const F64 B = JackelOneOverSqrtTwoPi * exp((-0.5 * (H * H + T * T))) * Expansion;
	return fabs(Max(B, 0.0));
}

internal F64
JackelBlackCallCdf(F64 X, F64 S)
{
	const F64 H = X / S, T = 0.5 * S, BMax = exp(0.5 * X),
			  B = JackelNormalCdf(H + T) * BMax - JackelNormalCdf(H - T) / BMax;
	return fabs(Max(B, 0.0));
}

internal F64
JackelNormalisedBlackCall(F64 X, F64 S)
{
	if (X > 0)
	{
		return JackelNormalisedIntrinsic(X, 1) + JackelNormalisedBlackCall(-X, S);
	}
	if (S <= 0)
	{
		return JackelNormalisedIntrinsic(X, 1);
	}
	if (X < S * JackelAsymptoticExpansionAccuracyThreshold &&
		0.5 * S * S + X <
			S * (JackelSmallTExpansionThreshold + JackelAsymptoticExpansionAccuracyThreshold))
	{
		return JackelBlackCallAsymptotic(X / S, 0.5 * S);
	}
	if (0.5 * S < JackelSmallTExpansionThreshold)
	{
		return JackelBlackCallSmallT(X / S, 0.5 * S);
	}
	if (X + 0.5 * S * S > S * 0.85)
	{
		return JackelBlackCallCdf(X, S);
	}
	return JackelBlackCallErfcx(X / S, 0.5 * S);
}

internal F64
JackelNormalisedVega(F64 X, F64 S)
{
	const F64 AbsX = fabs(X);
	return (AbsX <= 0)
		? JackelOneOverSqrtTwoPi * exp(-0.125 * S * S)
		: ((S <= 0 || S <= AbsX * JackelSqrtF64Min)
			? 0
			: JackelOneOverSqrtTwoPi *
				exp(-0.5 * (JackelSquare(X / S) + JackelSquare(0.5 * S))));
}

internal void
JackelLowerMap(F64 X, F64 S, F64 *F, F64 *Fp, F64 *Fpp)
{
	const F64 AbsX = fabs(X), Z = JackelSqrtOneOverThree * AbsX / S, Y = Z * Z, S2 = S * S,
			  Cdf = JackelNormalCdf(-Z), Pdf = JackelNormalPdf(Z);
	*Fpp = JackelPiOverSix * Y / (S2 * S) * Cdf *
		(8 * JackelSqrtThree * S * AbsX + (3 * S2 * (S2 - 8) - 8 * X * X) * Cdf / Pdf) *
		exp(2 * Y + 0.25 * S2);
	if (JackelIsZero(S))
	{
		*Fp = 1;
		*F = 0;
	}
	else
	{
		const F64 Cdf2 = Cdf * Cdf;
		*Fp = JackelTwoPi * Y * Cdf2 * exp(Y + 0.125 * S * S);
		*F = JackelIsZero(X) ? 0 : JackelTwoPiOverSqrtTwentySeven * AbsX * (Cdf2 * Cdf);
	}
}

internal F64
JackelInverseLowerMap(const F64 X, const F64 F)
{
	return JackelIsZero(F)
		? 0
		: fabs(X / (JackelSqrtThree * JackelInverseNormalCdf(
			pow(F / (JackelTwoPiOverSqrtTwentySeven * fabs(X)), 1. / 3.))));
}

internal void
JackelUpperMap(const F64 X, const F64 S, F64 *F, F64 *Fp, F64 *Fpp)
{
	*F = JackelNormalCdf(-0.5 * S);
	if (JackelIsZero(X))
	{
		*Fp = -0.5;
		*Fpp = 0;
	}
	else
	{
		const F64 W = JackelSquare(X / S);
		*Fp = -0.5 * exp(0.5 * W);
		*Fpp = JackelSqrtPiOverTwo * exp(W + 0.125 * S * S) * W / S;
	}
}

internal F64
JackelInverseUpperMap(F64 F)
{
	return -2. * JackelInverseNormalCdf(F);
}

internal F64
JackelUncheckedNormalisedImpliedVolatility(F64 Beta, F64 X, F64 Q, I32 MaxIterations)
{
	if (Q * X > 0)
	{
		Beta = fabs(Max(Beta - JackelNormalisedIntrinsic(X, Q), 0.));
		Q = -Q;
	}
	if (Q < 0)
	{
		X = -X;
		Q = -Q;
	}
	if (Beta <= 0)
	{
		return 0;
	}
	const F64 BMax = exp(0.5 * X);
	if (Beta >= BMax)
	{
		return DBL_MAX;
	}
	I32 Iterations = 0, DirectionReversalCount = 0;
	F64 F = -DBL_MAX, S = -DBL_MAX, Ds = S, DsPrevious = 0, SLeft = DBL_MIN, SRight = DBL_MAX;
	const F64 SC = sqrt(fabs(2 * X)), BC = JackelNormalisedBlackCall(X, SC),
			  VC = JackelNormalisedVega(X, SC);
	if (Beta < BC)
	{
		const F64 SL = SC - BC / VC, BL = JackelNormalisedBlackCall(X, SL);
		if (Beta < BL)
		{
			F64 FLowerMapL, DFLowerMapLDBeta, D2FLowerMapLDBeta2;
			JackelLowerMap(X, SL, &FLowerMapL, &DFLowerMapLDBeta, &D2FLowerMapLDBeta2);
			const F64 RLl = JackelConvexCubicControlRight(
				0., BL, 0., FLowerMapL, 1., DFLowerMapLDBeta, D2FLowerMapLDBeta2, 1);
			F = JackelRationalCubic(Beta, 0., BL, 0., FLowerMapL, 1., DFLowerMapLDBeta, RLl);
			if (!(F > 0))
			{
				const F64 T = Beta / BL;
				F = (FLowerMapL * T + BL * (1 - T)) * T;
			}
			S = JackelInverseLowerMap(X, F);
			SRight = SL;
			for (; Iterations < MaxIterations && fabs(Ds) > DBL_EPSILON * S; ++Iterations)
			{
				if (Ds * DsPrevious < 0)
				{
					++DirectionReversalCount;
				}
				if (Iterations > 0 && (3 == DirectionReversalCount || !(S > SLeft && S < SRight)))
				{
					S = 0.5 * (SLeft + SRight);
					if (SRight - SLeft <= DBL_EPSILON * S)
					{
						break;
					}
					DirectionReversalCount = 0;
					Ds = 0;
				}
				DsPrevious = Ds;
				const F64 B = JackelNormalisedBlackCall(X, S), Bp = JackelNormalisedVega(X, S);
				if (B > Beta && S < SRight)
				{
					SRight = S;
				}
				else if (B < Beta && S > SLeft)
				{
					SLeft = S;
				}
				if (B <= 0 || Bp <= 0)
				{
					Ds = 0.5 * (SLeft + SRight) - S;
				}
				else
				{
					const F64 LnB = log(B), LnBeta = log(Beta), BpOverB = Bp / B, H = X / S,
							  BHalley = H * H / S - S / 4,
							  Newton = (LnBeta - LnB) * LnB / LnBeta / BpOverB,
							  Halley = BHalley - BpOverB * (1 + 2 / LnB);
					const F64 BHh3 = BHalley * BHalley - 3 * JackelSquare(H / S) - 0.25,
							  Hh3 = BHh3 +
						2 * JackelSquare(BpOverB) * (1 + 3 / LnB * (1 + 1 / LnB)) -
						3 * BHalley * BpOverB * (1 + 2 / LnB);
					Ds = Newton * JackelHouseholderFactor(Newton, Halley, Hh3);
				}
				S += Ds = Max(-0.5 * S, Ds);
			}
			return S;
		}
		else
		{
			const F64 VL = JackelNormalisedVega(X, SL),
					  RLm = JackelConvexCubicControlRight(BL, BC, SL, SC, 1 / VL, 1 / VC, 0.0, 0);
			S = JackelRationalCubic(Beta, BL, BC, SL, SC, 1 / VL, 1 / VC, RLm);
			SLeft = SL;
			SRight = SC;
		}
	}
	else
	{
		const F64 SH = VC > DBL_MIN ? SC + (BMax - BC) / VC : SC,
				  BH = JackelNormalisedBlackCall(X, SH);
		if (Beta <= BH)
		{
			const F64 VH = JackelNormalisedVega(X, SH),
					  RHm = JackelConvexCubicControlLeft(BC, BH, SC, SH, 1 / VC, 1 / VH, 0.0, 0);
			S = JackelRationalCubic(Beta, BC, BH, SC, SH, 1 / VC, 1 / VH, RHm);
			SLeft = SC;
			SRight = SH;
		}
		else
		{
			F64 FUpperMapH, DFUpperMapHDBeta, D2FUpperMapHDBeta2;
			JackelUpperMap(X, SH, &FUpperMapH, &DFUpperMapHDBeta, &D2FUpperMapHDBeta2);
			if (D2FUpperMapHDBeta2 > -JackelSqrtF64Max && D2FUpperMapHDBeta2 < JackelSqrtF64Max)
			{
				const F64 RHh = JackelConvexCubicControlLeft(
					BH, BMax, FUpperMapH, 0., DFUpperMapHDBeta, -0.5, D2FUpperMapHDBeta2, 1);
				F = JackelRationalCubic(
					Beta, BH, BMax, FUpperMapH, 0., DFUpperMapHDBeta, -0.5, RHh);
			}
			if (F <= 0)
			{
				const F64 H = BMax - BH, T = (Beta - BH) / H;
				F = (FUpperMapH * (1 - T) + 0.5 * H * T) * (1 - T);
			}
			S = JackelInverseUpperMap(F);
			SLeft = SH;
			if (Beta > 0.5 * BMax)
			{
				for (; Iterations < MaxIterations && fabs(Ds) > DBL_EPSILON * S; ++Iterations)
				{
					if (Ds * DsPrevious < 0)
					{
						++DirectionReversalCount;
					}
					if (Iterations > 0 &&
						(3 == DirectionReversalCount || !(S > SLeft && S < SRight)))
					{
						S = 0.5 * (SLeft + SRight);
						if (SRight - SLeft <= DBL_EPSILON * S)
						{
							break;
						}
						DirectionReversalCount = 0;
						Ds = 0;
					}
					DsPrevious = Ds;
					const F64 B = JackelNormalisedBlackCall(X, S), Bp = JackelNormalisedVega(X, S);
					if (B > Beta && S < SRight)
					{
						SRight = S;
					}
					else if (B < Beta && S > SLeft)
					{
						SLeft = S;
					}
					if (B >= BMax || Bp <= DBL_MIN)
					{
						Ds = 0.5 * (SLeft + SRight) - S;
					}
					else
					{
						const F64 BMaxMinusB = BMax - B, G = log((BMax - Beta) / BMaxMinusB),
								  Gp = Bp / BMaxMinusB;
						const F64 BHalley = JackelSquare(X / S) / S - S / 4,
								  BHh3 = BHalley * BHalley - 3 * JackelSquare(X / (S * S)) - 0.25;
						const F64 Newton = -G / Gp, Halley = BHalley + Gp,
								  Hh3 = BHh3 + Gp * (2 * Gp + 3 * BHalley);
						Ds = Newton * JackelHouseholderFactor(Newton, Halley, Hh3);
					}
					S += Ds = Max(-0.5 * S, Ds);
				}
				return S;
			}
		}
	}
	for (; Iterations < MaxIterations && fabs(Ds) > DBL_EPSILON * S; ++Iterations)
	{
		if (Ds * DsPrevious < 0)
		{
			++DirectionReversalCount;
		}
		if (Iterations > 0 && (3 == DirectionReversalCount || !(S > SLeft && S < SRight)))
		{
			S = 0.5 * (SLeft + SRight);
			if (SRight - SLeft <= DBL_EPSILON * S)
			{
				break;
			}
			DirectionReversalCount = 0;
			Ds = 0;
		}
		DsPrevious = Ds;
		const F64 B = JackelNormalisedBlackCall(X, S), Bp = JackelNormalisedVega(X, S);
		if (B > Beta && S < SRight)
		{
			SRight = S;
		}
		else if (B < Beta && S > SLeft)
		{
			SLeft = S;
		}
		const F64 Newton = (Beta - B) / Bp, Halley = JackelSquare(X / S) / S - S / 4,
				  Hh3 = Halley * Halley - 3 * JackelSquare(X / (S * S)) - 0.25;
		S += Ds = Max(-0.5 * S, Newton * JackelHouseholderFactor(Newton, Halley, Hh3));
	}
	return S;
}

internal F64
JackelImpliedVolatility(F64 Price, F64 Forward, F64 Strike, F64 Expiry, option_type Type)
{
	F64 Q = (F64)Type;
	const F64 Intrinsic = fabs(Max(Q < 0 ? Strike - Forward : Forward - Strike, 0.0));
	if (Price < Intrinsic)
	{
		return -DBL_MAX;
	}
	const F64 MaxPrice = Q < 0 ? Strike : Forward;
	if (Price >= MaxPrice)
	{
		return DBL_MAX;
	}
	const F64 X = log(Forward / Strike);
	if (Q * X > 0)
	{
		Price = fabs(Max(Price - Intrinsic, 0.0));
		Q = -Q;
	}
	return JackelUncheckedNormalisedImpliedVolatility(
		Price / (sqrt(Forward) * sqrt(Strike)), X, Q, 2) / sqrt(Expiry);
}

internal void
JackelImpliedVolVec(const F64 *Prices, const F64 *Strikes, const option_type *Types,
	F64 Forward, F64 Expiry, F64 *Out, U64 N)
{
	for (U64 I = 0; I < N; I++)
	{
		Out[I] = JackelImpliedVolatility(Prices[I], Forward, Strikes[I], Expiry, Types[I]);
	}
}

