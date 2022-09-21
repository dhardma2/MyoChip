%MATLAB code to define peaks and troughs from waves in data from
%optogenetic stimulation of myofibers.
%Created by David Hardman, Usher Institute, University of Edinburgh
clear;
prompt = {'Input Filename:','Time Frame Length (s):','Pixel Width (um):','No. Rows;', 'Local max/min separation (Frames):','Output Filename(1st Wave):','Output Filename (mean values):'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'OPT025_03_Nucleus2_CTR.txt','1','1','98','10','OPT025_03_Nucleus2_CTR_Results_W1.txt','OPT025_03_Nucleus2_CTR_Results_mean_vals.txt'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

filename=char(answer(1));
Frame_Length=str2double(answer(2));
Pixel_width=str2double(answer(3));
Rows=str2double(answer(4));
t=str2double(answer(5));
t=t*10;
output_wave_1=char(answer(6));
output_mean=char(answer(7));


% Import data from text file
[FrameNr, Velocitymagnitudepx]=import_data(filename);
%Create spline
[curve, goodness, output]=fit(FrameNr(1:Rows),Velocitymagnitudepx(1:Rows),'smoothingspline','SmoothingParam',1);
x=1:(Rows*10);
for i=1:(Rows*10)
    j=i/10;
    y(i)=curve(j);
end
figure;plot(x,y)
%locate local maxima and minima, separated by t frames
LMax2 = islocalmax(y,'MinSeparation',t);
LMin2 = islocalmin(y,'MinSeparation',t);
%Sort maxima and minima to find apex and nadir of each wave
[Maxima_x, Maxima_y, Minima_x, Minima_y]=SortMaxMin(x,y,LMax2,LMin2);
xlim([0 length(x)])
hold on
scatter(Maxima_x, Maxima_y,'r*')
scatter(Minima_x, Minima_y,'b*')

function [FrameNr, Velocitymagnitudepx]=import_data(filename)
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
Velocitymagnitudepx = tbl.Velocitymagnitudepxframemean_valuemean4;
% Clear temporary variables
clear opts tbl
end

function [Maxima_x, Maxima_y, Minima_x, Minima_y]=SortMaxMin(x,y,LMax2,LMin2)
Raw_Maxima_x=x(LMax2);
Raw_Maxima_y=y(LMax2);
Raw_Minima_x=x(LMin2);
Raw_Minima_y=y(LMin2);

N_maxima=length(Raw_Maxima_x);
ii=0;
jj=0;
min_no=1;
if Raw_Minima_x(1)<Raw_Maxima_x(1)
    min_no=min_no+1;
end
max_no=1;
while max_no<=N_maxima-1
    k=0;
    for j=1:length(Raw_Minima_x)
        if Raw_Minima_x(j) > Raw_Maxima_x(max_no) && Raw_Minima_x(j) < Raw_Maxima_x(max_no+1)
            k=k+1;
            Temp_minima_x(k)=Raw_Minima_x(j);
            Temp_minima_y(k)=Raw_Minima_y(j);
        end
    end
    if k==0
        LMinTemp = islocalmin(y(Raw_Maxima_x(max_no):Raw_Maxima_x(max_no+1)));
        temp_xvals= x(LMinTemp)+Raw_Maxima_x(max_no);
        temp_yvals= y(temp_xvals);
        ii=ii+1;
        if Raw_Maxima_y(max_no)-temp_yvals(1)<0.05 %
            if Raw_Maxima_y(max_no)>Raw_Maxima_y(max_no+1)
                Maxima_y(ii)=Raw_Maxima_y(max_no);
                Maxima_x(ii)=Raw_Maxima_x(max_no);
                max_no=max_no+2;
            else
                Maxima_y(ii)=Raw_Maxima_y(max_no+1);
                Maxima_x(ii)=Raw_Maxima_x(max_no+1);
                max_no=max_no+2;
            end
        else
            Maxima_y(ii)=Raw_Maxima_y(max_no);
            Maxima_x(ii)=Raw_Maxima_x(max_no);
            max_no=max_no+1;
            jj=jj+1;
            Minima_y(jj)=temp_yvals(1);
            Minima_x(jj)=temp_xvals(1);
        end
    elseif k==1
        ii=ii+1;
        Maxima_y(ii)=Raw_Maxima_y(max_no);
        Maxima_x(ii)=Raw_Maxima_x(max_no);
        max_no=max_no+1;
        jj=jj+1;
        Minima_y(jj)=Raw_Minima_y(min_no);
        Minima_x(jj)=Raw_Minima_x(min_no);
        min_no=min_no+1;
    elseif k>1
        ii=ii+1;
        Maxima_y(ii)=Raw_Maxima_y(max_no);
        Maxima_x(ii)=Raw_Maxima_x(max_no);
        max_no=max_no+1;
        
        for L=1:k
            if Temp_minima_y(L)==min(Temp_minima_y)
                jj=jj+1;
                Minima_y(jj)=Temp_minima_y(L);
                Minima_x(jj)=Temp_minima_y(L);
                min_no=min_no+k;
            end
        end
    end
    clear Temp_minima_y Temp_minima_x
end
end
