//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <exception>
#include <functional>
#include <mutex>
#include <vector>

#include "../common/parallel.h"

//processInParallel that carries the first exception thrown by a task back
//to the calling thread (an exception escaping a std::thread would terminate)
template <class T>
void processInParallelRethrow(const std::vector<T>& tasks,
							  std::function<void(const T&)> updateFun,
							  size_t maxThreads, bool progressBar)
{
	std::exception_ptr error = nullptr;
	std::mutex errorMutex;
	std::function<void(const T&)> wrapped =
		[&updateFun, &error, &errorMutex] (const T& task)
	{
		{
			std::lock_guard<std::mutex> lock(errorMutex);
			if (error) return;
		}
		try
		{
			updateFun(task);
		}
		catch (...)
		{
			std::lock_guard<std::mutex> lock(errorMutex);
			if (!error) error = std::current_exception();
		}
	};
	processInParallel(tasks, wrapped, maxThreads, progressBar);
	if (error) std::rethrow_exception(error);
}
