
#include <deque>

typedef int dataType;

#define FS 360	// Sampling frequency.

	struct threshold
	{
		dataType i1 = 0;
		dataType i2 = 0;
		dataType spk = 0;
		dataType npk = 0;

		dataType thresh(const dataType signal, const dataType noise) const
				{ return noise + 0.25 * (signal - noise); }

		void set(dataType i)
		{
			i1 = i;
			i2 = 0.5 * i1;
		}

		void half()
		{
			i1 = 0.5 * i1;
			//i2 = 0.5 * i1; //???
		}

		void updateNoise(const dataType peak)
		{
			npk = 0.125 * peak + 0.875 * npk;
			set(thresh(spk, npk));
		}

		void updateSignal(const dataType peak)
		{
			spk = 0.125 * peak + 0.875 * spk;
			set(thresh(spk, npk));
		}

		void updateSignalSearchBack(const dataType peak)
		{
			spk = 0.25 * peak + 0.75 * spk;
			set(thresh(spk, npk));
		}
	};


	// DC filter

	std::deque<dataType> DCFilter(const std::deque<dataType> & signal)
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

	// Moving-Window Integration
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

	dataType slope(size_t i, std::deque<dataType> & squared) const
	{
		dataType currentSlope = 0;
		for (size_t j = i - 10; j <= i; j++)
			if (squared[j] > currentSlope) currentSlope = squared[j];
		return currentSlope;
	}

	struct panTompkins
	{
		size_t window;
		size_t rrmin;
		size_t rrmax;
		
		threshold threshold_i;
		threshold threshold_f;

		int rr1[8];
		int rr2[8];
		int rravg1;
		int rravg2;
		int rrlow = 0;
		int rrhigh = 0;
		int rrmiss = 0;

		bool regular = true;

	public:

	explicit panTompkins(size_t fs) :
		window(0.15*fs),
		rrmin(0.2*fs),
		rrmax(0.36*fs),
		
	{
	}

	std::deque<bool> detect(const std::deque<dataType> & signal);

	private:
		void updateRR(int index)
		{
			rravg1 = 0;
			for (size_t i = 0; i < 7; ++i)
			{
				rr1[i] = rr1[i+1];
				rravg1 += rr1[i];
			}
			rr1[7] = index; //add
			rravg1 += rr1[7];
			rravg1 *= 0.125; // 1/8
	
			if ((rr1[7] >= rrlow) && (rr1[7] <= rrhigh))
			{
				rravg2 = 0;
				for (size_t i = 0; i < 7; ++i)
				{
					rr2[i] = rr2[i+1];
					rravg2 += rr2[i];
				}
				rr2[7] = rr1[7]; //add
				rravg2 += rr2[7];
				rravg2 *= 0.125; // 1/8
	
				rrlow = 0.92 * rravg2;
				rrhigh = 1.16 * rravg2;
				rrmiss = 1.66 * rravg2;
			}

			bool prevRegular = regular;
			if (rravg1 == rravg2)
			{
				regular = true;
			} else {
				regular = false;
				if (prevRegular)
				{
					threshold_i.half();
					threshold_f.half();
				}
			}
		}

	};

std::deque<bool> panTompkins::detect(const std::deque<dataType> & signal)
{
	//TODO check for size, if more/less than XXX, throw exception

	// DC was not proposed on the original paper.
	// It is not necessary and can be removed if your sensor or database has no DC noise.
	std::deque<dataType> dcblock = DCFilter(signal);

	//1) Bandpass filter = LP + HP filter, filters work only for 200 Hz sampling rate
	std::deque<dataType> lowpass = lowPassFilter(dcblock); 
	std::deque<dataType> highpass = highPassFilter(lowpass);
	//normalize

	//TODO for another fs
	std::deque<dataType> bandpass = highpass;

	//2) Differentiator
	std::deque<dataType> derivative = derivativeFilter(bandpass);
	//normalize

	//3) Squaring
	std::deque<dataType> squared = squaredFilter(derivative);
	//4) Moving window integration
	std::deque<dataType> integral = MWI(squared, window);

	for (size_t i = 0; i < 8; ++i)
	{
		rr1[i] = 0;
		rr2[i] = 0;
	}

	long unsigned int i;

	long unsigned int lastQRS = 0;
	dataType lastSlope = 0;

	regular = true;

	std::deque<bool> rpeaks;

	for (long unsigned int current = 0; current < signal.size(); ++current)
	{
		bool qrs = false;

		long unsigned int sample = current+1; //???

		if ((integral[current] >= threshold_i.i1) &&
		    (bandpass[current] >= threshold_f.i1))
		{
			if (sample - lastQRS > rrmin)
			{
				dataType currentSlope = slope(current, squared);
				if (sample - lastQRS <= rrmax)
				{
					if (currentSlope <= lastSlope/2)
					{
						qrs = false;
					} else
					{
						lastSlope = currentSlope;
						threshold_i.updateSignal(integral[current]);
						threshold_f.updateSignal(bandpass[current]);
						qrs = true;
					}
				} else
				{
					lastSlope = currentSlope;
					threshold_i.updateSignal(integral[current]);
					threshold_f.updateSignal(bandpass[current]);
					qrs = true;
				}
			} else
			{
				threshold_i.updateNoise(integral[current]);
				threshold_f.updateNoise(bandpass[current]);
				qrs = false;
				rpeaks[current] = qrs;
				//if (sample > DELAY + BUFFSIZE) output(rpeaks[0]);
				continue;
			}
		}

		if (qrs)
		{
			updateRR(sample - lastQRS);
			lastQRS = sample;
		} else
		{
			//back search
			if ((sample - lastQRS > (long unsigned int)rrmiss) &&
			    (sample - lastQRS > rrmin))
			{
				//do search
				for (i = current - (sample - lastQRS) + rrmin;
					i < (long unsigned int)current;
					++i)
				{
					if ( (integral[i] > threshold_i.i2) &&
					     (bandpass[i] > threshold_f.i2))
					{
						dataType currentSlope = slope(i, squared);
						if (
							(currentSlope < lastSlope/2) &&
							(i + sample) < lastQRS + 0.36*lastQRS)
						{
							qrs = false;
						} else
						{
							lastSlope = currentSlope;
							threshold_i.updateSignalBackSearch(integral[i]);
							threshold_f.updateSignalBackSearch(bandpass[i]);
							qrs = true;

							updateRR(sample - (current - i) - lastQRS);
							lastQRS = sample - (current - i);
							break;
						}
					}
				} //for

				if (qrs)
				{
					rpeaks[current] = false;
					rpeaks[i] = true;
					//if (sample > DELAY + BUFFSIZE) output(rpeaks[0]);
					continue;
				}
			} //if

			if (!qrs)
			{
				if ((integral[current] >= threshold_i.i1) ||
				    (bandpass[current] >= threshold_f.i1))
				{
					threshold_i.updateNoise(integral[current]);
					threshold_f.updateNoise(bandpass[current]);
				}
			}
		} // !qrs

		rpeaks[current] = qrs;
		//if (sample > DELAY + BUFFSIZE) output(rpeaks[0]);
	}

	//TODO Q,S peaks
	return rpeaks;
}

