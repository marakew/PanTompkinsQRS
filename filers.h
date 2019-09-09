#ifndef _fillers_h_
#define _fillers_h_
#endif
	template<typename dataType>
	void normalize(std::deque<dataType> & values)
	{
		dataType max_value = values[0];
		for (size_t index = 1; index < values.size(); ++index)
			if (values[index] > max_value) max_value = values[index];
	
		for (size_t index = 0; index < values.size(); ++index)
			values[i] /= max_value;
	}

	// DC filter
	template<typename dataType>
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
	template<typename dataType>
	std::deque<dataType> lowPassFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = signal[index];
			if (index >= 1) value += 2 * result[index - 1];
			if (index >= 2) value -= result[index - 2];
			if (index >= 6) value -= 2 * signal[index - 6];
			if (index >= 12) value += signal[index - 12];
			result.push_back(value);
		}
		return result;
	}

	// High Pass filter
	// Implemented as proposed by the original paper.
	// y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]
	template<typename dataType>
	std::deque<dataType> highPassFilter(const std::deque<dataType> & signal)
	{
		std::deque<dataType> result;
		for (size_t index = 0; index < signal.size(); ++index)
		{
			dataType value = -signal[index];
			if (index >= 1) value -= result[index - 1];
			if (index >= 16) value += 32 * signal[index - 16];
			if (index >= 32) value += signal[index - 32];
			result.push_back(value);
		}
		return result;
	}


	struct BandFilter
	{
		double b0, b1, b2,
			a1, a2;
		double x0, x1, x2,
			y1, y2;

		reset()
		{
			x0 = 0; x1 = 0; x2 = 0; y1 = 0; y2 = 0;
		}

		BandFilter(double cutoff, bool hp)
		{
			const double B = tan(cutoff * M_PI);
			const double BB = B * B;
			const double S = 1.0 + M_SQRT2 * B + BB;

			if (hp) {
				b0 = 1.0 / S;
				b1 = -2.0 * b0;
			} else {
				b0 = BB / S;
				b1 = 2.0 * b0;
			}
			b2 = b0;
			a1 = 2.0 * (BB - 1.0) / S;
			a2 = (1.0 - M_SQRT2 * B + BB) / S;			
		}

		void filter(std::deque<dataType> & signal)
		{
			for (size_t i = 0; i < signal.size(); ++i)
			{
				x2 = x1;
				x1 = x0;
				x0 = signal[i];
				signal[i] = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
				y2 = y1;
				y1 = signal[i];
			}
		}
	};

	template<typename dataType>
	std::deque<dataType> bandPassFilter(const std::deque<dataType> & signal, double lowcut, double highcut, double rate)
	{
		std::deque<dataType> result = signal;
		BandFilter lowPass(lowcut/rate, false);
		lowPass.filter(result);
		BandFilter highPass(highcut/rate, true);
		highPass.filter(result);
		return result;
	}

	// Derivative filter
	// Implemented as proposed by the original paper.
	// y(nT) = (1/8T)[-x(nT - 2T) - 2x(nT - T) + 2x(nT + T) + x(nT + 2T)]
	template<typename dataType>
	std::deque<dataType> derivativeFilter(const std::deque<dataType> & signal, size_t fs = 1)
	{
		double T = 1.0/fs; //???
		std::deque<dataType> result;
		for (size_t index = 2; index < signal.size() - 2; ++index)
		{
			dataType value = -signal[index - 2] - 2 * signal[index - 1] +
					2 * signal[index + 1] + signal[index + 2];
			value /= (8.0 * T);
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
	template<typename dataType>
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
	template<typename dataType>
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
