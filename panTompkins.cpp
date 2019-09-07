
#include <memory>
#include <stdexcept>
#include <deque>
#include "threshold.h"
#include "filters.h"

	typedef int dataType;

	const std::deque<size_t> detectQpeaks(const std::deque<dataType> & bandpass, const std::deque<size_t> & rpeaks);
	const std::deque<size_t> detectSpeaks(const std::deque<dataType> & bandpass, const std::deque<size_t> & rpeaks);

	struct panTompkins
	{
	private:
		size_t samplefrequency;
		size_t window;
		size_t rrmin;
		size_t rrmax;
		size_t slopemax;
		
		threshold<dataType> threshold_i;
		threshold<dataType> threshold_f;

		int rr1[8];
		int rr2[8];
		int rravg1;
		int rravg2;
		int rrlow = 0;
		int rrhigh = 0;
		int rrmiss = 0;

		bool regular = true;
	public:
		std::deque<size_t> qpeaks;
		std::deque<size_t> rpeaks;
		std::deque<size_t> speaks;

	explicit panTompkins(size_t fs = 360) :
		samplefrequency(fs),
		window(0.15*fs), //150 ms
		rrmin(0.2*fs), //200 ms
		rrmax(0.36*fs), //360 ms
		slopemax(0.75*fs), //10 ???
		
	{
	}

		void detectPeaks(const std::deque<dataType> & signal);

	private:
		dataType findMaxSlope(size_t i, std::deque<dataType> & values)
		{
			assert(i>=slopemax);
			dataType max_value = 0;
			for (size_t j = i - slopemax; j <= i; ++j)
				if (values[j] > max_value) max_value = values[j];
			return max_value;
		}

		void detectRpeaks(const std::deque<dataType> & bandpass);

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

void panTompkins::detectPeaks(const std::deque<dataType> & signal)
{
	if (signal.size() < samplefrequency*2)
		throw std::length_error("input signal too short");

	// DC was not proposed on the original paper.
	// It is not necessary and can be removed if your sensor or database has no DC noise.
	std::deque<dataType> dcblock = dcFilter(signal);

	//1) Bandpass filter = LP + HP filter, filters work only for 200 Hz sampling rate
	std::deque<dataType> lowpass = lowPassFilter(dcblock); 
	std::deque<dataType> highpass = highPassFilter(lowpass);
	//??? normalize(highpass);

	//TODO for another fs
	std::deque<dataType> bandpass = highpass;

	//qpeaks.clear();
	//rpeaks.clear();
	//speaks.clear();

	//TODO calculate sum delay by filters and flush peaks from begin

	detectRpeaks(bandpass);

	qpeaks = detectQpeaks(bandpass, rpeaks); //left Min
	speaks = detectSpeaks(bandpass, rpeaks); //right Min
}

const std::deque<size_t> detectQpeaks(const std::deque<dataType> & bandpass, const std::deque<size_t> & peaks)
{
	std::deque<size_t> result;
	for (size_t i = 0; i < peaks.size(); ++i)
	{
		size_t n = peaks[i];
		if (n-1<0) break;
		while(bandpass[n]>bandpass[n-1])
		{
			--n;
			if (n<0) break;
		}
		result.push_back(n);
	}
	return result;
}

const std::deque<size_t> detectSpeaks(const std::deque<dataType> & bandpass, const std::deque<size_t> & peaks)
{
	std::deque<size_t> result;
	for (size_t i = 0; i < peaks.size(); ++i)
	{
		size_t n = peaks[i];
		if (n+1>=bandpass.size()) break;
		while(bandpass[n]>bandpass[n+1])
		{
			++n;
			if (n>=bandpass.size()) break;
		}
		result.push_back(n);
	}
	return result;
}

void panTompkins::detectRpeaks(const std::deque<dataType> & bandpass)
{
	//2) Differentiator
	std::deque<dataType> derivative = derivativeFilter(bandpass);
	//??? normalize(derivative);

	//3) Squaring
	std::deque<dataType> squared = squaredFilter(derivative);
	//4) Moving window integration
	std::deque<dataType> integral = MWI(squared, window);

	for (size_t i = 0; i < 8; ++i)
	{
		rr1[i] = 0;
		rr2[i] = 0;
	}

	size_t lastQRS = 0;
	dataType lastSlope = 0;

	regular = true;

	rpeaks.clear();

	std::deque<bool> peaks;

	for (size_t current = 0; current < bandpass.size(); ++current)
	{
		bool qrs = false;

		size_t sample = current+1; //???

		if ((integral[current] >= threshold_i.i1) &&
		    (bandpass[current] >= threshold_f.i1))
		{
			if (sample - lastQRS > rrmin)
			{
				dataType currentSlope = findMaxSlope(current, squared); //integral
				if (sample - lastQRS <= rrmax)
				{
					if (currentSlope <= lastSlope/2)
					{
						//T-wave found
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
				peaks[current] = qrs;
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
			if ((sample - lastQRS > (size_t)rrmiss) &&
			    (sample - lastQRS > rrmin))
			{
				size_t i;
				//do search
				for (i = current - (sample - lastQRS) + rrmin;
					i < (size_t)current;
					++i)
				{
					if ( (integral[i] > threshold_i.i2) &&
					     (bandpass[i] > threshold_f.i2))
					{
						dataType currentSlope = findMaxSlope(i, squared); //integral
						if (
							(currentSlope < lastSlope/2) &&
							(i + sample) < lastQRS + 0.36*lastQRS)
						{
							//T-wave found
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
					peaks[current] = false;
					peaks[i] = true;
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

		peaks[current] = qrs;
		//if (sample > DELAY + BUFFSIZE) output(rpeaks[0]);
	}

	for (size_t i = 0; i < peaks.size(); ++i)
	{
		if (peaks[i])
			rpeaks.push_back(i);
	}
}

