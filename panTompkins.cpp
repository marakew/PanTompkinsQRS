
#include <memory>
#include <stdexcept>
#include <deque>
#include "threshold.h"
#include "filters.h"

	typedef int dataType;

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
		std::deque<bool> qpeaks;
		std::deque<bool> rpeaks;
		std::deque<bool> speaks;

	explicit panTompkins(size_t fs = 360) :
		samplefrequency(fs),
		window(0.15*fs),
		rrmin(0.2*fs),
		rrmax(0.36*fs),
		slopemax(0.75*fs),
		
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

		void detectQpeaks(const std::deque<dataType> & bandpass);
		void detectRpeaks(const std::deque<dataType> & bandpass);
		void detectSpeaks(const std::deque<dataType> & bandpass);

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

	qpeaks.clear();
	rpeaks.clear();
	speaks.clear();

	//TODO calculate sum delay by filters and flush rpeaks from begin by false values

	detectRpeaks(bandpass);

	detectQpeaks(bandpass); //left Min
	detectSpeaks(bandpass);	//right Min
}

void panTompkins::detectQpeaks(const std::deque<dataType> & bandpass)
{
	//TODO
}

void panTompkins::detectSpeaks(const std::deque<dataType> & bandpass)
{
	//TODO
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

	for (size_t current = 0; current < bandpass.size(); ++current)
	{
		bool qrs = false;

		size_t sample = current+1; //???

		if ((integral[current] >= threshold_i.i1) &&
		    (bandpass[current] >= threshold_f.i1))
		{
			if (sample - lastQRS > rrmin)
			{
				dataType currentSlope = findMaxSlope(current, squared);
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
						dataType currentSlope = findMaxSlope(i, squared);
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
}

