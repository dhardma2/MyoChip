clear;
prompt = {'Pixel width (microns):','Day:','Position:','Enhancement Actin:','Enhancement Nuclei:','image length:','Sample number:','Record striations?:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.3405','1','1','5','15','500','10','yes'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

D=str2double(answer(2));
%pixel size
ps=str2double(answer(1));
pos=str2double(answer(3));
EFA=str2double(answer(4));
EFN=str2double(answer(5));
L1=str2double(answer(6));
L2=round(L1/ps,0);
SampSize=str2double(answer(7));
StriCheck = 0;
if strcmpi(answer(8), 'yes')
    StriCheck = 1;
end
filename1=sprintf('pos%d-D%d-Actin.tif',pos,D);
Actin=imread(filename1)*EFA;

filename2=sprintf('pos%d-D%d-DNA.tif',pos,D);
Nuclei=imread(filename2)*EFN;

[MeanNucDist,SDNucDist,Coeff_var,stri,totl]=ManDist(Actin,Nuclei,L2,SampSize,StriCheck);
MeanNucDist = MeanNucDist * ps;
SDNucDist = SDNucDist * ps;
GlobalMeanNucDist=mean(MeanNucDist(MeanNucDist>0 & SDNucDist >0));
GlobalSDNucDist= mean(SDNucDist(MeanNucDist>0 & SDNucDist >0));
GlobalCoeff_var = mean(Coeff_var(Coeff_var>0));
GlobalTot_length = mean(totl);

    
clear Actin Nuclei

function [MeanNucDist,SDNucDist,Coeff_var,stri,totl]=ManDist(Actin,Nuclei,len,SampSize,StriCheck)
j = 0;
C=imfuse(Nuclei(1:len,1:len), Actin(1:len,1:len),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);  
figure;imshow(C)
title('Click on consecutive nuclei in myotube nearest * marker. Press space') 
xsamp=rand(1,SampSize)*len;
ysamp=rand(1,SampSize)*len;
hold on
scatter(xsamp,ysamp,'w*');
x = 0;
y = 0;
button = 1;
i = 0;
k = 0;
 xt1=0;
 yt1=0;
while k < SampSize
    [xg,yg,button] = ginput(1);
    
    if button==1
        i = i+1;
        xt1(i) = xg; yt1(i)= yg;
        hold on
        plot(xt1(i),yt1(i),'w.','MarkerSize',20)
        drawnow
    end
    if button ==32
            k=k+1;
            j=j+1;
            pt(:,1)=xt1; pt(:,2)=yt1;
            d1=0;
            for n=1:length(xt1)-1
                d1(n)=norm(pt(n+1,:)-pt(n,:));
            end
                MeanNucDist(j)=mean(d1);
                totl(j)=sum(d1);
                if sum(d1) >0
                    SDNucDist(j)=std(d1,1);
                    Coeff_var(j)= SDNucDist(j)/MeanNucDist(j);
                else
                    SDNucDist(j)=0;
                    Coeff_var(j)= 0;
                end
                if StriCheck == 1          
                    list = {'Full Striations','Partial Striations','No striations','Skip',};
                    [indx] = listdlg('ListString',list);
                    if indx ==1
                     stri(j)=1;
                 elseif indx ==2
                        stri(j)=2;
                    elseif indx ==3
                        stri(j)=0;
                    else
                        stri(j)=3;
                    end
                end
                i = 0;
                xt1=0;%clear xt1
                yt1=0;%clear yt1
                clear pt
   end
    
                
end

end
