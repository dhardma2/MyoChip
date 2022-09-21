%MATLAB code for calculating statistics from optogenetic stimulation of 
%myofiber data.
%Written by David Hardman, Usher Institute, University of Edinburgh

Maxima_x=Maxima_x/10;
Minima_x=Minima_x/10;

%Stats for first twitch
Accel_Disp_Wave_1_Raw=Maxima_y(1);
Accel_Disp_Wave_1=Accel_Disp_Wave_1_Raw*Pixel_width;

Accel_Time_1_Raw=(Maxima_x(1));
Accel_Time_1=Accel_Time_1_Raw*Frame_Length;

Accel_Rate_Wave_1_Raw=Accel_Disp_Wave_1_Raw/Accel_Time_1_Raw;
Accel_Rate_Wave_1=Accel_Rate_Wave_1_Raw*Pixel_width/Frame_Length;

Decel_Disp_Wave_1_Raw=Maxima_y(1)-Minima_y(1);
Decel_Disp_Wave_1=Decel_Disp_Wave_1_Raw*Pixel_width;

Decel_Time_1_Raw=Minima_x(1)-Maxima_x(1);
Decel_Time_1=Decel_Time_1_Raw*Frame_Length;

Decel_Rate_Wave_1_Raw=Accel_Disp_Wave_1_Raw/Decel_Time_1_Raw;
Decel_Rate_Wave_1=Decel_Rate_Wave_1_Raw*Pixel_width/Frame_Length;

Wavelength_Wave_1_Raw=Minima_x(1);
Wavelength_Wave_1=Wavelength_Wave_1_Raw*Frame_Length;

Frequency_Wave_1=1/Wavelength_Wave_1;
Tfirst=table(Accel_Disp_Wave_1,Accel_Time_1,Accel_Rate_Wave_1,Decel_Disp_Wave_1,Decel_Time_1,Decel_Rate_Wave_1,Wavelength_Wave_1);

%average stats for subsequent twitches
for wave=1:length(Maxima_x)-1
    wave1=wave+1;
    Accel_Disp(wave)=Maxima_y(wave1)-Minima_y(wave);
    Accel_Time(wave)=Maxima_x(wave1)-Minima_x(wave);
    Accel_Rate(wave)=Accel_Disp(wave)/Accel_Time(wave);
    
    Decel_Disp(wave)=Maxima_y(wave1)-Minima_y(wave1);
    Decel_Time(wave)=Minima_x(wave1)-Maxima_x(wave1);
    Decel_Rate(wave)=Decel_Disp(wave)/(Minima_x(wave1)-Maxima_x(wave1));
    
    Wavelength(wave)=Minima_x(wave1)-Minima_x(wave);
    Frequency(wave)=1/(Wavelength(wave)*Frame_Length);
end

Mean_Accel_Disp=mean(Accel_Disp);
Mean_Accel_Disp=Mean_Accel_Disp*Pixel_width;

SD_Accel_Disp=std(Accel_Disp);
SD_Accel_Disp=SD_Accel_Disp*Pixel_width;

Mean_Accel_Time=mean(Accel_Time);
Mean_Accel_Time=Mean_Accel_Time*Frame_Length;

SD_Accel_Time=std(Accel_Time);
SD_Accel_Time=SD_Accel_Time*Frame_Length;

Mean_Decel_Disp=mean(Decel_Disp);
Mean_Decel_Disp=Mean_Decel_Disp*Pixel_width;

SD_Decel_Disp=std(Decel_Disp);
SD_Decel_Disp=SD_Decel_Disp*Pixel_width;

Mean_Decel_Time=mean(Decel_Time);
Mean_Decel_Time=Mean_Decel_Time*Frame_Length;

SD_Decel_Time=std(Decel_Time);
SD_Decel_Time=SD_Decel_Time*Frame_Length;

Mean_Accel_Rate=mean(Accel_Rate);
Mean_Accel_Rate=Mean_Accel_Rate*Pixel_width/Frame_Length;

SD_Accel_Rate=std(Accel_Rate);
SD_Accel_Rate=SD_Accel_Rate*Pixel_width/Frame_Length;

Mean_Decel_Rate=mean(Decel_Rate);
Mean_Decel_Rate=Mean_Decel_Rate*Pixel_width/Frame_Length;

SD_Decel_Rate=std(Decel_Rate);
SD_Decel_Rate=SD_Decel_Rate*Pixel_width/Frame_Length;

Mean_Wavelength=mean(Wavelength);
Mean_Wavelength=Mean_Wavelength*Frame_Length;

SD_Wavelength=std(Wavelength);
SD_Wavelength=SD_Wavelength*Frame_Length;

Mean_Frequency=mean(Frequency);

SD_Frequency=std(Frequency);

Tmean=table(Mean_Accel_Disp,SD_Accel_Disp,Mean_Accel_Time,SD_Accel_Time,Mean_Decel_Disp,SD_Decel_Disp,Mean_Decel_Time,SD_Decel_Time,Mean_Accel_Rate,SD_Accel_Rate,Mean_Decel_Rate,SD_Decel_Rate,Mean_Wavelength,SD_Wavelength,Mean_Frequency,SD_Frequency);

writetable(Tfirst,output_wave_1)

writetable(Tmean,output_mean)
