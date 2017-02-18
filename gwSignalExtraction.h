#include "gwDataTypes.h"

class Extractor
{
	public:
		Extractor();

		//Setup memeber functions
		void setSignal(Signal* input);
		void setTemplates(std::vector<Template>* templates);

		//Overloaded function to perform fft on a 
		//set of tenplates or signals
		void fft(std::vector<Template>* templatesFFT);
		void fft(Signal* signalFFT);
		void fftComplex(std::vector<Template>* templatesFFT);
		void fftComplex(Signal* signalFFT);
		
		void fftInverse(std::vector<Signal>* signalsFFTI);
		void fftInverseComplex(std::vector<Signal>* signalsFFTI);
		
		//Performs time/frequency domain convolutions of the respective signal 
		//against all templates in the respective template bank
		void tConvolution(std::vector<Signal>* output);
		void fConvolution(std::vector<Signal>* output);

	private:
		Signal* mSignalT;
		Signal* mSignalF;

		std::vector<Template>* mTemplatesT;
		std::vector<Template>* mTemplatesF;
};
