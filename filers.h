#ifndef _fillers_h_
#define _fillers_h_
#endif
	void normalize(std::deque<dataType> & values)
	{
		dataType max_value = values[0];
		for (size_t index = 1; index < values.size(); ++index)
			if (values[index] > max_value) max_value = values[index];
	
		for (size_t index = 0; index < values.size(); ++index)
			values[i] /= max_value;
	}

	// DC filter

	std::deque<dataType> dcFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = 0;
			if (index >= 1) value = signal[index] - signal[index - 1] + 0.995 * result[index - 1];
			result.push_back(value);
		}
		return result;
	}

	// Low Pass filter
	// Implemented as proposed by the original paper.
	// y(nT) = 2y(nT - T) - y(nT - 2T) + x(nT) - 2x(nT - 6T) + x(nT - 12T)

	std::deque<dataType> lowPassFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = signal[index];
			if (index >= 1) value += 2 * result[index - 1];
			if (index >= 2) value -= result[index - 2];
			if (index >= 6) value -= 2 * result[index - 6];
			if (index >= 12) value += result[index - 12];
			result.push_back(value);
		}
		return result;
	}

	// High Pass filter
	// Implemented as proposed by the original paper.
	// y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]

	std::deque<dataType> highPassFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = -signal[index];
			if (index >= 1) value -= result[index - 1];
			if (index >= 16) value += 32 * result[index - 16];
			if (index >= 32) value += result[index - 32];
			result.push_back(value);
		}
		return result;
	}


	// Derivative filter
	// Implemented as proposed by the original paper.
	// y(nT) = (1/8T)[-x(nT - 2T) - 2x(nT - T) + 2x(nT + T) + x(nT + 2T)]

	std::deque<dataType> derivativeFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 2; index < signal.size() - 2; ++index)
		{
			dataType value = -signal[index - 2] - 2 * signal[index - 1] +
					2 * signal[index + 1] + signal[index + 2];
			value /= 8.0;
			result.push_back(value);
		}
		dataType value;
		value = result.front();
		result.push_front(value);
		result.push_front(value);
		value = result.back();
		result.push_back(value);
		result.push_back(value);
		return result;
	}

	// Squared filter
	// Implemented as proposed by the original paper.
	// y(nT) = [x(nT)]^2.

	std::deque<dataType> squaredFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = signal[index]*signal[index];
			result.push_back(value);
		}
		return result;
	}

	// Moving-Window Integration, delay N/2 ???
	// Implemented as proposed by the original paper.
	// y(nT) = (1/N)[x(nT - (N - 1)T) + x(nT - (N - 2)T) + ... x(nT)]
	// N, in samples, must be defined so that the window is ~150ms.

	std::deque<dataType> MWI(const std::deque<dataType> & signal, size_t N)
	{
		std::deque<dataType> result;
		dataType value = 0;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			int first = index - (N - 1); //???
			value += signal[index] / N;
			if (first > 0) value -= signal[first - 1] / N;
			result.push_back(value);
		}
		return result;
	}

#endif
