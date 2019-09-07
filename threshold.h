#ifndef _threshold_h_
#define _threshold_h_

	inline dataType thresh(const dataType signal, const dataType noise) 
		{ return noise + 0.25 * (signal - noise); }

	template<typename dataType>
	struct threshold
	{
		dataType i1 = 0;
		dataType i2 = 0;
		dataType spk = 0;
		dataType npk = 0;

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

		void update()
		{
			set(thresh(spk, npk));
		}

		void updateNoise(const dataType peak)
		{
			npk = 0.125 * peak + 0.875 * npk;
			update();
		}

		void updateSignal(const dataType peak)
		{
			spk = 0.125 * peak + 0.875 * spk;
			update();
		}

		void updateSignalSearchBack(const dataType peak)
		{
			spk = 0.25 * peak + 0.75 * spk;
			update();
		}
	};

#endif
