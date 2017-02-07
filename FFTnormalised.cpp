#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <cmath>


std::vector<double> FFTNormalised(std::vector<double>);//declare functions
std::vector<double> FrequencyConvertor(std::vector<double>);
void outputWriter(std::vector<double>, std::vector<double>);

int main(){
	
	std::ifstream myfile;//open a file reader
	myfile.open ("Q1_a_time_domain_dataset.dat");//open the file
	std::vector<double> fTime;//declare vector arrays
	std::vector<double> fAmplitude;
	double a,b;//these hold the read in data so that the vectors can be pushed back
	//int i =0;	//set counter
while(!myfile.eof()){//while more text in the file
	
	myfile>>a>>b;//read in values of the row
	fTime.push_back(a);//put data into arrays 
	fAmplitude.push_back(b);
	//std::cout<<fTime[i]<<"Ampl  "<<fAmplitude[i]<<std::endl;
	//i++;
} 

	
	std::vector<double> fourierAmplitude=FFTNormalised(fAmplitude);	

	std::vector<double> frequency = FrequencyConvertor(fTime);

	outputWriter(frequency,fourierAmplitude);
	return 0;
}


void outputWriter(std::vector<double>Frequency,std::vector<double> fourierAmplitude){//This function writes out the two vectors two a file which can be plotted in MATLAB
	
	std::ofstream newfile;
	newfile.open("Outputforgraphing.txt");
	const int m=fourierAmplitude.size()/2;

for(int j=0;j<m;j++){
	
	newfile<<Frequency[j]<<"       "<<fourierAmplitude[j]<<std::endl;
	//std::cout<<fTime[i]<<"Ampl  "<<fAmplitude[i]<<std::endl;
	//i++;
}
}
std::vector<double> FFTNormalised(std::vector<double> Amplitude){//This function produces a normalised Fourier transform of the amplitudes
	
	int n=Amplitude.size();//declare how much data is being transformed
	std::vector<double> AmplitudeHolder=Amplitude;
	gsl_fft_real_wavetable * real;//declare the use of the real waveforms
	//gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work;//declare an area for the function to use
	
	work = gsl_fft_real_workspace_alloc (n);//assign a work of area
	real = gsl_fft_real_wavetable_alloc (n);
	
	gsl_fft_real_transform (&AmplitudeHolder[0], 1, n, real, work);//tranform the data
                               
	gsl_fft_real_wavetable_free (real);//cleear the wavetable
	
	const int m=Amplitude.size()/2;
	std::vector<double> fourierAmplitude;
	
	for(int l=0;l<m;l++){//this normalises the Amplitude vector
		fourierAmplitude.push_back(AmplitudeHolder[l]/(double)m);
	}
	
	return fourierAmplitude;
}
std::vector<double> FrequencyConvertor(std::vector<double> Time){//converts a time array into a frequency domain
	int n=Time.size();	
	double samplefreq=n/Time[n-1];//
	const int m=n/2;
	std::vector<double> freq;
	
	
	for(int k=0;k<m;k++){
		freq.push_back((double)k*samplefreq/(2*(double)n));//create a vector that is a match for the fourierAmplitude vector
	}
	return freq;
}





