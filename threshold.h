#ifndef _threshold_h_
#define _threshold_h_

	template<typename dataType>
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
			//i2 = 0.5 * i1; //??? imp - no, algo - yes
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

#endif
