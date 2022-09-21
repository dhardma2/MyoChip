clear;
%Matlab code for performing Fast Fourier Transform on PIV data, isolating
%the dominant frequency and providing statistics for comparison.
% Written by David Hardman, Usher Institute, University of Edinburgh

prompt = {'Input Filename:','Time Frame Length (s):','Pixel Width (um):','No. Rows;','Output Filename:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'OPT027_09_Nucleus2_TRAINED.txt','0.001','1','98','OPT027_09_Nucleus2_TRAINED_FFT.txt'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

filename=char(answer(1));
Frame_Length=str2double(answer(2));
Pixel_width=str2double(answer(3));
Rows=str2double(answer(4));
output_FFT=char(answer(5));
%Determine sampling frequency
Fs=1/Frame_Length;
% Import data from text file
[FrameNr, X]=import_data(filename,Rows);
%Convert to actual distance and time
X=X*Pixel_width;
t=FrameNr*Frame_Length;
%Fast Fourier Transform of data
Y=fft(X);

L=length(FrameNr);

%Calculate magnitudes from FFT
P2 = abs(Y/L);

P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
%figure;plot(P1)
f = Fs*(0:(L/2))/L;
figure;plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%Find frequency with maximum magnitude (dominant frequency)
P1_2=P1;
P1_2(1)=0;
[F_Max,I] = max(P1_2);
Y_phase=angle(Y);

%frequency
freq_Max_Mag= f(I);
%wavelenth
lambda_max_mag= 1/freq_Max_Mag;
%Length of acceleration phase
Accel_phase_max_mag=lambda_max_mag/2;
%Displacement during acceleration phase
disp_max_mag=2*F_Max;
%Phase angle
phase_max_mag=Y_phase(I);

%Determine data for other frequencies with local maxima
P1_maxima=islocalmax(P1_2);
P1_maxima=P1_maxima.*P1;
P1filt=P1_maxima>=F_Max/2;

freq_all=f(P1filt);
lambda_all=1./freq_all;
disp_all=2*P1_2(P1filt);
disp_all=disp_all';

phase_all=Y_phase(P1filt);

%Plot original data against dominant frequency wave
figure;plot(t,X)
S1=P1(1)+F_Max*sin(2*pi*freq_Max_Mag*(t)+phase_max_mag);

hold on
plot(t,S1)
title('PIV waveform and dominant frequency')
xlabel('time (s)')
ylabel('Displacement (\mum)')

%Record and write data to file
T_FFT=table(freq_Max_Mag,lambda_max_mag,Accel_phase_max_mag,disp_max_mag,phase_max_mag);
writetable(T_FFT,output_FFT)

function [FrameNr, Velocitymagnitudepx]=import_data(filename,Rows)
% Import data from text file
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["FrameNr4", "Velocitymagnitudepxframemean_valuemean4"];
opts.VariableTypes = ["double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
tbl = readtable(filename, opts);
% Convert to output type
FrameNr = tbl.FrameNr4;
FrameNr=FrameNr(1:Rows);
Velocitymagnitudepx = tbl.Velocitymagnitudepxframemean_valuemean4;
Velocitymagnitudepx =Velocitymagnitudepx(1:Rows);
% Clear temporary variables
clear opts tbl
end