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
std::vector<double> matchedFilter(std::vector<double>, std::vector<double>);
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


	std::vector<double> fAreas=matchedFilter(fTime,  fAmplitude);

	
	//outputWriter(runfreq,fAreas);
	
	return 0;

}

std::vector<double> matchedFilter(std::vector<double> Time, std::vector<double> Amplitude){//this function return the amplitude of a matched filter of a simple sin wave in the time domain
	
	double timeIncrement=Time[Time.size()-1]/(Time.size());//calculate and declare necessary constants
	double sampleFrequency=1/timeIncrement;
	double area;
	const int siz2=(int)(sampleFrequency/2);
	std::vector<double> runfreq;
	std::vector<double> Areas;
	const int size=Time.size();
	std::vector<double> Y;
	std::vector<double> innerproduct;
	
	for(int i =0; i<siz2;i++){//check every frequency

		runfreq.push_back(i);//set run frequency to interger value
		area=0.0;//set area for this run
	
		for(int j=0; j<size;j++){//calculate the intergrand () the data multiplied by the template)

			Y.push_back(sin(runfreq[i]*Time[j]*2*pi));
			innerproduct.push_back(Y[j]*Amplitude[j]);
		}
	
		for(int k =0; k<size-1; k++){
		
			area+=timeIncrement*(.5*(innerproduct[k]+innerproduct[k+1]));//calcualte the intergral
		
		}
	
		
		Areas.push_back(area);//add this area to the array
		Y.clear();//clear vectors for next run
		innerproduct.clear();//
	
	}
	
	outputWriter(runfreq,Areas);
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



