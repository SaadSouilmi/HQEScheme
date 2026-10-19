#ifndef RNG_H
#define RNG_H

#if !defined(__SIZEOF_INT128__)
#error "PCG64 requires compiler support for unsigned 128-bit integers."
#endif
typedef unsigned __int128 U128;
struct random_series { U128 State, Stream; };

// Tables are built once, safely across threads. RandomSeed also calls this.
internal void RandomInit(void);

// PCG XSL-RR 128/64; Seed and Stream use the low 64 bits of PCG's seed inputs.
internal void RandomSeed(random_series *R, U64 Seed, U64 Stream);
internal U64 RandomU64(random_series *R);
internal U32 RandomU32(random_series *R);
internal F64 RandomUniform(random_series *R); // (0,1), 52 random bits; HQE/log(U).
internal F64 RandomUniformClosedOpen(random_series *R); // [0,1), 53 bits; Ziggurat/log1p(-U).
internal F64 RandomNormal(random_series *R);
internal F64 RandomExponential(random_series *R);


#endif // RNG_H

