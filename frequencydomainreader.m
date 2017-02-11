
%%
fAllData=csvread('Output1.csv');
%import 'Q1_a_time_domain_dataset.t.*'

fSize=size(fAllData);

for n=1:(fSize(2))
   fTime(n)=fAllData(1,n);
   fAmplitude(n)=fAllData(2,n); 
   fAmplitude2(n)=fAllData(3,n);
   fAmplitude3(n)=fAllData(4,n);
end
subplot(3,2,1)
plot(fTime,fAmplitude,'-')
subplot(3,2,2)
plot(fTime,fAmplitude2,'-')
subplot(3,2,3)
plot(fTime,fAmplitude3,'-')
%%
fAllData=csvread('Output2.csv');
%import 'Q1_a_time_domain_dataset.t.*'

fSize=size(fAllData);

for n=1:(fSize(2))
  % fTime2(n)=fAllData(1,n);
   fAmplitude4(n)=fAllData(2,n); 
   fAmplitude5(n)=fAllData(3,n);
   fAmplitude6(n)=fAllData(1,n);
end
subplot(3,2,4)
plot(fTime,fAmplitude4,'-')
subplot(3,2,5)
plot(fTime,fAmplitude5,'-')
subplot(3,2,6)
plot(fTime,fAmplitude6,'-')

%fId = fopen('fre.csv','r')







