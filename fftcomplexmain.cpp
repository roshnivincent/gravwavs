
#include "gwSignalExtraction.h"
#include "gwReadWrite.h"

int main(){

	std::string file1="NonAdaptiveTemplates.dat";
	std::vector<Template> Templates;
	delimiter f=tab;
	
	loadTemplates(file1,&Templates,f);
		
	std::string file2="NonAdaptiveSignal.dat";
	
	std::vector<Signal> signal;
		
	loadSignals(file2,&signal,f);
	
	Extractor Bob;
	
	Bob.setSignal(&signal[0]);
	Bob.setTemplates(&Templates);
	
	Bob.fftComplex(&signal[0]);
	//Bob.fftComplex(&Templates);
	
	
	
	


	return 0;
}