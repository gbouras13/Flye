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
struct PortableRandState
{
	std::mutex mtx;
	uint32_t r[34];
	uint64_t n = 0;
};

inline PortableRandState& portableRandState()
{
	static PortableRandState state;
	return state;
}

inline int portableRand()
{
	PortableRandState& s = portableRandState();
	std::lock_guard<std::mutex> lock(s.mtx);
	if (s.n == 0)
	{
		int32_t word = 1;	//glibc default seed
		s.r[0] = word;
		for (int i = 1; i < 31; ++i)
		{
			//word = (16807 * word) % 2147483647 without overflow
			int32_t hi = word / 127773;
			int32_t lo = word % 127773;
			word = 16807 * lo - 2836 * hi;
			if (word < 0) word += 2147483647;
			s.r[i] = word;
		}
		for (int i = 31; i < 34; ++i) s.r[i] = s.r[i - 31];
		s.n = 34;
		for (; s.n < 344; ++s.n) s.r[s.n % 34] = s.r[(s.n - 31) % 34] + s.r[(s.n - 3) % 34];
	}
	uint32_t val = s.r[(s.n - 31) % 34] + s.r[(s.n - 3) % 34];
	s.r[s.n % 34] = val;
	++s.n;
	return (int)(val >> 1);
}

//restarts the sequence, as in a freshly started process
//(used when several modules run in one process)
inline void portableRandReset()
{
	PortableRandState& s = portableRandState();
	std::lock_guard<std::mutex> lock(s.mtx);
	s.n = 0;
}
