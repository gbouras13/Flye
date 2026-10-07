//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <cstdint>
#include <mutex>

//Platform-independent replacement for rand(). The C library rand() differs
//between glibc, macOS libc and musl (e.g. WebAssembly builds), so the reads
//sampled for parameter estimation - and hence the assembly - depended on the
//platform. This reproduces the sequence of glibc's default-seeded rand(),
//so results match Flye built on Linux.
inline int portableRand()
{
	static std::mutex mtx;
	static uint32_t r[34];
	static uint64_t n = 0;

	std::lock_guard<std::mutex> lock(mtx);
	if (n == 0)
	{
		int32_t word = 1;	//glibc default seed
		r[0] = word;
		for (int i = 1; i < 31; ++i)
		{
			//word = (16807 * word) % 2147483647 without overflow
			int32_t hi = word / 127773;
			int32_t lo = word % 127773;
			word = 16807 * lo - 2836 * hi;
			if (word < 0) word += 2147483647;
			r[i] = word;
		}
		for (int i = 31; i < 34; ++i) r[i] = r[i - 31];
		n = 34;
		for (; n < 344; ++n) r[n % 34] = r[(n - 31) % 34] + r[(n - 3) % 34];
	}
	uint32_t val = r[(n - 31) % 34] + r[(n - 3) % 34];
	r[n % 34] = val;
	++n;
	return (int)(val >> 1);
}
