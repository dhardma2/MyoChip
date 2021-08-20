clearvars -except MCmean0
%input image variables
prompt = {'Pixel width (microns):','Day:','Position:','Enhancement factor:','Image length:','Dilation factor:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.3405','1','1','20','500','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

D=str2double(answer(2));
%pixel size
ps=str2double(answer(1));
pos=str2double(answer(3));
EF=str2double(answer(4));
L1=str2double(answer(5));
DilF=str2double(answer(6));
L2=round(L1/ps,0);
%Read images, binarize images and crop to length
filename1=sprintf('pos%d-D%d-Actin.tif',pos,D);
Actin=imread(filename1)*EF;
Actin=Actin(1:L2,1:L2);
MTBW = imbinarize(Actin,'global');
filename2=sprintf('pos%d-D%d-DNA.tif',pos,D);
Nuclei=imread(filename2)*EF;
Nuclei=Nuclei(1:L2,1:L2);
MCBW = imbinarize(Nuclei,'global');

%open small regions of the image to remove remaining myotubes/artifacts and dilate the remaining regions. 
BWactin2=bwareaopen(MTBW,2000);
seMT = strel('disk',DilF);
BWactinDil=imdilate(BWactin2,seMT);
MTNuc=MCBW & BWactinDil;
MCNuc=~MCBW | MTNuc;
MCNuc=~MCNuc;
statsMC = regionprops('table',MCNuc,'Centroid','area');
%calculate a mean myocyte area for estimating numbers of myocytes in clusters.
if D==1
    MCmean=mean(statsMC.Area);
    MCmean0=MCmean;
else
    emc = exist('MCmean0');
    if emc == 1
        MCmean=MCmean0;
    else
        MCmean=mean(statsMC.Area);
    end
end
%remove areas which are smaller than 3/4 of the mean area.
MCNuc=bwareaopen(MCNuc,ceil(3*MCmean/4));

%figure; imshow(imoverlay(MCNuc,MTNuc))
clear statsMC;
statsMC=regionprops('table',MCNuc,'Centroid','area');
statsMT=regionprops('table',MTNuc,'Centroid','area');
MTTotArea=sum(statsMT.Area);
MTEst=MTTotArea/MCmean;
x=statsMC.Centroid(:,1);
y=statsMC.Centroid(:,2);
xT=statsMT.Centroid(:,1);
yT=statsMT.Centroid(:,2);
C=imfuse(Nuclei(1:L2,1:L2), Actin(1:L2,1:L2),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);  
figure(1);imshow(C);

hold on
scatter(x,y,L2/5,'w.')
scatter(xT,yT,L2/5,'b.')
%figure(1);plot(xt1(:,1), yt1(:,1), 'r.', 'MarkerSize', 10)

%Function to manually label cells missed by segmentation 
[x1,y1,x2,y2]=ManLab(Actin);
MCTot=length(x)-length(x1)+length(x2);
MTTot=round(MTEst+length(x1)-length(x2),0);

%nuclei size correction based on mean
MTArea=statsMT.Area;
i=0;
j=0;
for a=1:length(MTArea)
    if MTArea(a)<= 2*MCmean
        i=i+1;
        MTSingNuc(i)=MTArea(a);
    else 
        j=j+1;
        MTMultiNuc(j)=MTArea(a);
    end
end

i=0;
for a=1:length(MTSingNuc)
    if MTSingNuc(a)>= 100
        i=i+1;
        MTSingNucLg(i)=MTSingNuc(a);
    end
end
MTMultiEst=sum(MTMultiNuc)/(2*mean(MTSingNucLg));
MTEst2=length(MTSingNuc)+ MTMultiEst+length(x1)-length(x2);
MTEst2=round(MTEst2,0);
FIndex=MTEst2/(MCTot+MTEst2);
clear MCBW MCNuc MTBW MTNuc Actin BWactin2 BWactinDil Nuclei C ...
    a answer D i j MTEst MTMultiEst MTMultiNuc MTTot x x1 x2 y1 y2


function [xt1,yt1,xt2,yt2]=ManLab(Iin)
title('Left click->mytube, Right click -> myocyte, spacebar-> finish')
[sizex]=size(Iin(:,1,:,1));
[sizey]=size(Iin(1,:,:,1));
for row = 1 : 500 : sizex(1)%3154
 line([1, sizey(2)], [row, row]);
end
for column = 1 : 500 : sizey(2)%3954
 line([column , column ], [1, sizex(1)]);
end
xt1 = 0;
yt1 = 0;
xt2 = 0;
yt2 = 0;
x = 0;
y = 0;
button = 1;
i = 0;
j = 0;
while button <=3
    [xg,yg,button] = ginput(1);
    
    if button==1
        j = j+1;
        xt1(j) = xg; yt1(j)= yg;
        hold on
        figure(1);plot(xt1(j),yt1(j),'b.','MarkerSize',20)
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
    
    if button==3
        i = i+1;
        xt2(i) = xg; yt2(i)= yg;
        hold on
        figure(1);plot(xt2(i),yt2(i),'w.','MarkerSize',20)
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
end

end