//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Helpers that reproduce Python semantics used by Flye's Python pipeline,
//so that the C++ port (flye-modules pipeline) gives identical results.

#pragma once

#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace pycompat
{

//CPython's random.Random (MT19937, seeded with init_by_array as in
//_randommodule.c), with getrandbits / _randbelow / shuffle as in random.py
class PyRandom
{
public:
	explicit PyRandom(uint32_t seed)
	{
		uint32_t key[1] = {seed};
		this->initByArray(key, 1);
	}

	uint32_t genrandUint32()
	{
		static const uint32_t MAG01[2] = {0x0U, 0x9908b0dfU};
		uint32_t y;
		if (_index >= N)
		{
			int kk;
			for (kk = 0; kk < N - M; ++kk)
			{
				y = (_mt[kk] & UPPER_MASK) | (_mt[kk + 1] & LOWER_MASK);
				_mt[kk] = _mt[kk + M] ^ (y >> 1) ^ MAG01[y & 0x1U];
			}
			for (; kk < N - 1; ++kk)
			{
				y = (_mt[kk] & UPPER_MASK) | (_mt[kk + 1] & LOWER_MASK);
				_mt[kk] = _mt[kk + (M - N)] ^ (y >> 1) ^ MAG01[y & 0x1U];
			}
			y = (_mt[N - 1] & UPPER_MASK) | (_mt[0] & LOWER_MASK);
			_mt[N - 1] = _mt[M - 1] ^ (y >> 1) ^ MAG01[y & 0x1U];
			_index = 0;
		}
		y = _mt[_index++];
		y ^= (y >> 11);
		y ^= (y << 7) & 0x9d2c5680U;
		y ^= (y << 15) & 0xefc60000U;
		y ^= (y >> 18);
		return y;
	}

	//random.getrandbits(k) for 1 <= k <= 32
	uint32_t getrandbits(int k)
	{
		if (k < 1 || k > 32) throw std::runtime_error("getrandbits: unsupported k");
		return this->genrandUint32() >> (32 - k);
	}

	//random._randbelow_with_getrandbits(n), n > 0
	uint64_t randbelow(uint64_t n)
	{
		int k = 0;
		for (uint64_t v = n; v; v >>= 1) ++k;	//n.bit_length()
		if (k > 32) throw std::runtime_error("randbelow: n too large");
		uint64_t r = this->getrandbits(k);
		while (r >= n) r = this->getrandbits(k);
		return r;
	}

	//random.shuffle(x)
	template <class T>
	void shuffle(std::vector<T>& x)
	{
		if (x.size() < 2) return;
		for (size_t i = x.size() - 1; i >= 1; --i)
		{
			size_t j = (size_t)this->randbelow(i + 1);
			std::swap(x[i], x[j]);
		}
	}

private:
	static const int N = 624;
	static const int M = 397;
	static const uint32_t UPPER_MASK = 0x80000000U;
	static const uint32_t LOWER_MASK = 0x7fffffffU;

	void initGenrand(uint32_t s)
	{
		_mt[0] = s;
		for (int i = 1; i < N; ++i)
		{
			_mt[i] = (1812433253U * (_mt[i - 1] ^ (_mt[i - 1] >> 30)) + (uint32_t)i);
		}
		_index = N;
	}

	void initByArray(const uint32_t* key, size_t keyLength)
	{
		this->initGenrand(19650218U);
		size_t i = 1;
		size_t j = 0;
		size_t k = ((size_t)N > keyLength ? (size_t)N : keyLength);
		for (; k; --k)
		{
			_mt[i] = (_mt[i] ^ ((_mt[i - 1] ^ (_mt[i - 1] >> 30)) * 1664525U))
					 + key[j] + (uint32_t)j;
			++i; ++j;
			if (i >= (size_t)N) { _mt[0] = _mt[N - 1]; i = 1; }
			if (j >= keyLength) j = 0;
		}
		for (k = N - 1; k; --k)
		{
			_mt[i] = (_mt[i] ^ ((_mt[i - 1] ^ (_mt[i - 1] >> 30)) * 1566083941U))
					 - (uint32_t)i;
			++i;
			if (i >= (size_t)N) { _mt[0] = _mt[N - 1]; i = 1; }
		}
		_mt[0] = 0x80000000U;
	}

	uint32_t _mt[N];
	int _index;
};

//Python's round(x) for a float: round half to even
inline double pyRound(double x)
{
	double r = std::floor(x);
	double diff = x - r;
	if (diff > 0.5) return r + 1;
	if (diff < 0.5) return r;
	return (std::fmod(r, 2.0) == 0.0) ? r : r + 1;
}

//flye.utils.utils.get_median: middle element, or mean of the two middle ones
template <class T>
double getMedian(std::vector<T> values)
{
	if (values.empty()) throw std::runtime_error("_get_median() arg is an empty sequence");
	std::sort(values.begin(), values.end());
	size_t n = values.size();
	if (n % 2 == 1) return (double)values[n / 2];
	return ((double)values[n / 2 - 1] + (double)values[n / 2]) / 2;
}

//Insertion-ordered dictionary (Python dict semantics: re-assignment keeps
//the original position of the key)
template <class V>
class OrderedDict
{
public:
	V& operator[](const std::string& key)
	{
		auto it = _index.find(key);
		if (it != _index.end()) return _items[it->second].second;
		_index[key] = _items.size();
		_items.emplace_back(key, V());
		return _items.back().second;
	}

	bool contains(const std::string& key) const {return _index.count(key) > 0;}

	const V& at(const std::string& key) const
	{
		auto it = _index.find(key);
		if (it == _index.end()) throw std::out_of_range("OrderedDict: no key " + key);
		return _items[it->second].second;
	}

	V& at(const std::string& key)
	{
		auto it = _index.find(key);
		if (it == _index.end()) throw std::out_of_range("OrderedDict: no key " + key);
		return _items[it->second].second;
	}

	size_t size() const {return _items.size();}
	bool empty() const {return _items.empty();}

	typedef typename std::vector<std::pair<std::string, V>>::iterator iterator;
	typedef typename std::vector<std::pair<std::string, V>>::const_iterator const_iterator;
	iterator begin() {return _items.begin();}
	iterator end() {return _items.end();}
	const_iterator begin() const {return _items.begin();}
	const_iterator end() const {return _items.end();}

	std::vector<std::string> sortedKeys() const
	{
		std::vector<std::string> keys;
		for (auto& kv : _items) keys.push_back(kv.first);
		std::sort(keys.begin(), keys.end());
		return keys;
	}

private:
	std::vector<std::pair<std::string, V>> _items;
	std::unordered_map<std::string, size_t> _index;
};

//splits on runs of ASCII whitespace, like Python's str.split()
inline std::vector<std::string> pySplit(const std::string& s)
{
	std::vector<std::string> out;
	size_t i = 0;
	auto isSpace = [](char c)
		{return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';};
	while (i < s.size())
	{
		while (i < s.size() && isSpace(s[i])) ++i;
		if (i >= s.size()) break;
		size_t j = i;
		while (j < s.size() && !isSpace(s[j])) ++j;
		out.push_back(s.substr(i, j - i));
		i = j;
	}
	return out;
}

//splits on a single separator, like Python's str.split(sep)
inline std::vector<std::string> pySplitSep(const std::string& s, char sep)
{
	std::vector<std::string> out;
	size_t start = 0;
	while (true)
	{
		size_t pos = s.find(sep, start);
		if (pos == std::string::npos)
		{
			out.push_back(s.substr(start));
			break;
		}
		out.push_back(s.substr(start, pos - start));
		start = pos + 1;
	}
	return out;
}

//Python's str.strip() for ASCII whitespace
inline std::string pyStrip(const std::string& s)
{
	auto isSpace = [](char c)
		{return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';};
	size_t b = 0, e = s.size();
	while (b < e && isSpace(s[b])) ++b;
	while (e > b && isSpace(s[e - 1])) --e;
	return s.substr(b, e - b);
}

}	//namespace pycompat
