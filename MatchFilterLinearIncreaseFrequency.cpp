#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <stdio.h>

#define pi (atan(1)*4)

std::vector<double> FFTNormalised(std::vector<double>);//declare functions
std::vector<double> FrequencyConvertor(std::vector<double>);
std::vector<double> matchedFilterLinearFrequency(std::vector<double>, std::vector<double>);
void outputWriter(std::vector<double>, std::vector<double>);

int main(){
	
	std::vector<double> fTime,fAmplitude,fFrequency;//declare vector arrays

	
	for(int i=0;i<300000;i++){
		fTime.push_back(i/10000.0);
		fFrequency.push_back(30.0+fTime[i]*0.52);
		fAmplitude.push_back(10*sin(2*pi*fFrequency[i]*fTime[i]));
	}
	


	std::vector<double> fAreas=matchedFilterLinearFrequency(fTime, fAmplitude);
	//outputWriter(fTime,fAmplitude);
	
	

	
	return 0;

}

std::vector<double> matchedFilterLinearFrequency(std::vector<double> Time, std::vector<double> Amplitude){//this function finds a linearly increasing frequency term and plots it in the "increase increment domain"
	
	double timeIncrement=Time[Time.size()-1]/(Time.size());//calculate and declare necessary constants
	double sampleFrequency=1/timeIncrement;
	double area;
	std::vector<double> increaseinc;
	const int siz2=(int)(sampleFrequency/2);
	std::vector<double> runfreq;
	std::vector<double> Areas;
	const int size=Time.size();
	std::vector<double> Y;
	std::vector<double> innerproduct;
	double startfreq=30.0;//this function requires a start frequency 
	double currentfrequency;
	
	for(int i =0; i<100;i++){//check every increment requires a update to stop hardcoding

		//startfreq.push_back(60);//set run frequency to interger value
		area=0.0;//set area for this run
		increaseinc.push_back(0.01*i);//calculate the increase incremnt for this run
		
		for(int j=0; j<size;j++){//calculate the intergrand ( the data multiplied by the template)
			currentfrequency=increaseinc[i]*Time[j]+startfreq;//calculate what the current frequency 
			Y.push_back(sin(currentfrequency*Time[j]*2*pi));//create this stage of the template
			innerproduct.push_back(Y[j]*Amplitude[j]);//calcualte the product of template and data
		}
	
		for(int k =0; k<size-1; k++){
		
			area+=timeIncrement*(.5*(innerproduct[k]+innerproduct[k+1]));//calcualte the intergral
		
		}
	
		
		Areas.push_back(area);//add this area to the array
		Y.clear();//clear vectors for next run
		innerproduct.clear();//
	
	}
	
	outputWriter(increaseinc,Areas);
	return Areas;//return the amplitudes
}



void outputWriter(std::vector<double>Frequency,std::vector<double> fourierAmplitude){//This function writes out the two vectors to a file which can be plotted in MATLAB
	
	std::ofstream newfile;
	newfile.open("Outputforgraphing.csv");
	const int m=fourierAmplitude.size();

for(int j=0;j<m;j++){
	
	newfile<<Frequency[j]<<","<<fourierAmplitude[j]<<std::endl;
	//std::cout<<fTime[i]<<"Ampl  "<<fAmplitude[i]<<std::endl;
	//i++;
}
}



