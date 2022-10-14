function FeatureVector = FeatureCalc(RawData,config)
%% this function calculates the features of EEG signal used in Fan et al. 2018 Eur J Neurosci.
%   7 temporal features: absolute histogram of 7 equally sized bins
%   7 sepctral features: repative spectrum within 7 bands (0.5-1.5,1.5-2.5,2.5-4.5,4.5-8.5,8.5-16.5,16.5-32.6,32.5-64.5)
% inpupt:
%          - RawData: one-channel EEG signal
%		  - config: configuration of algorithm
%		      config.Fs: sample frequency;

LengthData = length(RawData);
N = LengthData;
Nfft = 2500;
Fs = config.Fs;
NormData = (RawData - mean(RawData))/std(RawData);
ProcessedData = NormData;

band=[0.5 1.5 2.5 4.5 8.5 16.5 32.5 64.5]; % frequency bands
%%  Symetric Histogram proposed by Jonanthan, number of bins are determined by number of frequency bands
h_Data = hist(abs(NormData),length(band)-1)/N;

%% Relative power in subbands using fft
NumUniquePts = ceil((Nfft+1)/2); 
f_resolution = Fs/Nfft;
DataFFT_all= abs(fft(NormData(1:N).*hamming(N),Nfft));
DataFFT_all = DataFFT_all(1:NumUniquePts);
DataFFT_all(1) = DataFFT_all(1)/N;
DataFFT_all(2:end)=DataFFT_all(2:end)*2/N;
Ptot = sum(DataFFT_all.^2);
for j=1:length(band)-1;
    band_init=ceil(band(j)/f_resolution);
    band_end=ceil(band(j+1)/f_resolution);
    Pmodel(1,j)=sum(DataFFT_all(band_init:band_end).^2)/Ptot;
end   
%% Relative power in subbands using fft
FeatureVector = [Pmodel,h_Data];
end