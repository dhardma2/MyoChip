%Code for obtaining myonuclei and myocyte counts and density from
%pericentrin and DAPI nuclei staining.
%Includes estimate of cell width and BTX staining

    
clearvars -except MCmean0
%input image variables
prompt = {'Pixel width (microns):','Day:','Start Position:','Enhancement factor:','image length:','Dilation factor:','% Pericentrin threshold:'};
dlgtitle = 'Input';
dims = [1 35];
%Pericentrin threshold will change with fluorescence intensity
definput = {'0.227','2','32','20','500','3','0.05'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

D=str2double(answer(2));
%pixel size
ps=str2double(answer(1));
pos=str2double(answer(3));
EF=str2double(answer(4));
L1=str2double(answer(5));
DilF=str2double(answer(6));
Pericentrin_Thresh=str2double(answer(7));
L2=round(L1/ps,0);
count_no=0;
for p=pos:pos+9
count_no=count_no+1;
%Read images, binarize images and crop to length
filename1=sprintf('2 Hz-Pos%d-MTDNA.tif',p);
filename2=sprintf('2 Hz-Pos%d-DNABW.tif',p);
filename3=sprintf('2 Hz-Pos%d-DNA.tif',p);
filename4=sprintf('2 Hz-Pos%d-MTDNABW.tif',p);
filename5=sprintf('2 Hz-Pos%d-ActinBW.tif',p);
filename6=sprintf('2 Hz-Pos%d-BTX.tif',p);

%Define myotube regions from actin staining and remove small regions
Myotube=imread(filename5);
Myotube=imbinarize(Myotube);
Myotube=Myotube(1:L2,1:L2);
Myotube=bwareaopen(Myotube,2000);
%Read Pericentrin stained regions for observation
Pericentrin=imread(filename1)*5*EF;
Pericentrin=Pericentrin(1:L2,1:L2);
%Read DAPI stained nuclei regions
NucBW=imread(filename2);
NucBW=NucBW(1:L2,1:L2);
NucBW=imbinarize(NucBW);
seMB = strel('disk',DilF);
%Dilate Nuclei
NucDil=imdilate(NucBW,seMB);
%Read unbinarized DAPI nuclei staining for observation
Nuclei=imread(filename3)*EF;
Nuclei=Nuclei(1:L2,1:L2);
%Read pericentrin stained regions
Pericentrin_Nuclei=imread(filename4);
Pericentrin_Nuclei=Pericentrin_Nuclei(1:L2,1:L2);
Pericentrin_Nuclei=imbinarize(Pericentrin_Nuclei);

Pericentrin_Nuclei=imdilate(Pericentrin_Nuclei,seMB);
%Calculate overlapping regions of Pericentrin and DAPI (pericentrin inside
%nuclei)
Pericentrin_Nuclei=Pericentrin_Nuclei & NucDil;
%figure;imshow(MBNuc);
%Centroid/area stats for DAPI stained nuclei regions
statsNuc = regionprops('table',NucBW,'Centroid','area');
   
%Define average nuclei size to use as a threshold for excluding smaller
%regions.
emc = exist('MCmean0');
if emc == 1
   MCmean=MCmean0;
else
   MCmean=mean(statsNuc.Area);
end
%remove areas which are smaller than 1/2 of the mean nucleus area.
Nuc=bwareaopen(NucBW,ceil(1*MCmean/2));
%figure; imshow(imoverlay(Nuc,NucBW))
clear statsMC;
%statsNuc=regionprops('table',Nuc,'Centroid','area');
%Choose whether to dilate the nuclei image again? 
%NucDil=imdilate(Nuc,seMB);
%NucDil=Nuc;
%Search through nuclei regions and count the number of pericentrin positive
%pixels
StatsNucPx=regionprops('table',Nuc,'PixelList','Centroid','area');
pxcount=0;
test=0;
for px1=1:length(StatsNucPx.PixelList)
    Pixels1=cell2mat(StatsNucPx.PixelList(px1));
    for px2=1:length(Pixels1)
        if logical(Pericentrin_Nuclei(Pixels1(px2,2),Pixels1(px2,1)))
            test(Pixels1(px2,2),Pixels1(px2,1))=1;
            pxcount=pxcount+1;
        end
    end
    %Calculate proportion of cell with pericentrin staining
    Proportion_pericentrin(px1)=pxcount/length(Pixels1);
    pxcount=0;
end       


%nuclei size correction based on mean
NucArea=StatsNucPx.Area;
i=0;
j=0;
MTMultiNuc=0;
for a=1:length(NucArea)
    if NucArea(a)<= 2*MCmean
        i=i+1;
        MTSingNuc(i)=NucArea(a);
    else 
        j=j+1;
        MTMultiNuc(j)=NucArea(a);
    end
end

i=0;
for a=1:length(MTSingNuc)
    if MTSingNuc(a)>= 100
        i=i+1;
        MTSingNucLg(i)=MTSingNuc(a);
    end
end

NucX=StatsNucPx.Centroid(:,1);
NucY=StatsNucPx.Centroid(:,2);
%Show overlay of DAPI and Pericentrin
C=imfuse(Nuclei(1:L2,1:L2), Pericentrin(1:L2,1:L2),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);  
figure(1);imshow(C);

hold on

x=0;

for k=1:px1
    if Proportion_pericentrin(k)<=Pericentrin_Thresh
        scatter(NucX(k),NucY(k),L2/10,'b.')
        x=x+1;
    else
        scatter(NucX(k),NucY(k),L2/10,'w.')
    end
end
%Function to manually add/remove myoblasts missed by segmentation 
title('Left click->add myoblast, Right click -> Remove myoblast, spacebar-> finish')
AddMBx=0;
AddMBy=0;
RemMBx=0;
RemMBy=0;
%[AddMBx,AddMBy,RemMBx,RemMBy]=ManLab(Pericentrin);
if sum(AddMBx)>0
    AddMB=length(AddMBx);
else
    AddMB=0;
end
if sum(RemMBx)>0
    RemMB=length(RemMBx);
else
    RemMB=0;
end
%Total myoblasts
MCTot=x-RemMB+AddMB;
Myoblast_estimate(count_no)=MCTot;
%Myoblasts to subtract from total nuclei count
MCSub=x-RemMB;
AddMTx=0;
AddMTy=0;
RemMTx=0;
RemMTy=0;
title('Left click->Remove myotube, Right click -> Add myotube, spacebar-> finish')
%[RemMTx,RemMTy,AddMTx,AddMTy]=ManLab(Pericentrin);
if sum(AddMTx)>0
    AddMT=length(AddMTx);
else
    AddMT=0;
end
if sum(RemMTx)>0
    RemMT=length(RemMTx);
else
    RemMT=0;
end

%Calculate Myotube density as a fraction of the sample area

MTsize=size(Myotube);
s1=MTsize(1);
s2=MTsize(2);

[MTDens,MTD]=MTDensity(Myotube,s1,s2);
Proportional_Area(count_no)=MTD;

%remove holes from myotube image
Myotube=imdilate(Myotube,seMB);
Filled=imfill(Myotube,'holes');
holes = Filled & ~Myotube;
bigholes = bwareaopen(holes, 200);
smallholes = holes & ~bigholes;
I = Myotube | smallholes;
%Calculate a reference Myotube width by dividing density by total length.
%I=Myotube(:,:);
%figure;imshow(I);
BWskel = bwskel(I,'MinBranchLength',50);
%figure;imshow(BWskel);
TotLength=sum(sum(BWskel));
ReferenceWidth=MTD/TotLength;
ReferenceWidth=ReferenceWidth*ps;
Reference_Width(count_no)=ReferenceWidth;

%Calculate mean intensity of BTX per myonucleus
%Find myonuclei by subtracting nuclei from myotubes
Myonuclei=I & NucBW;
Myonuclei=bwareaopen(Myonuclei,0.25*MCmean);
%figure;imshow(Myonuclei);
%Remove groups of myonuclei
LgMyonuclei=bwareaopen(Myonuclei,2*MCmean);
Myonuclei= Myonuclei & ~LgMyonuclei;
%figure;imshow(Myonuclei);
BTX=imread(filename6);
BTX=BTX(1:L2,1:L2);
%convert to 8 bit
%BTX=uint8(BTX);
StatsBTXNucPx=regionprops(Myonuclei,'PixelIdxList');
for k=1:numel(StatsBTXNucPx)
idx=StatsBTXNucPx(k).PixelIdxList;
MeanInt(k)=mean(BTX(idx));
end
matrixname=sprintf('Mean_BTX_Intensity-2HZ-Pos%d',p);
writematrix(MeanInt,matrixname);
%Estimate of area of multiple nuclei by dividing total area by mean of
%larger individual nuclei
MTMultiEst=sum(MTMultiNuc)/(mean(MTSingNucLg));
%Sum of single nuclei, estimated multi-nuclei and manually added myonuclei 
%minus manually removed nuclei and nuclei labelled as myocytes.
MTEst2=length(MTSingNuc)+AddMT - RemMT + MTMultiEst-MCSub;
MTEst2=round(MTEst2,0);
MyonucleiEstimate(count_no)=MTEst2;
%Fusion index
FIndex(count_no)=MTEst2/(MCTot+MTEst2);
clear C MBNuc Nuc NucBW NucDil Pericentrin Nuclei Pericentrin_Nuclei I Myotube BWskel test...
    holes bigholes smallholes filled BTX Myonuclei LgMyonuclei MeanInt
end
function [MTDens,MTD]=MTDensity(Myotubes,s1,s2)
MTD=sum(sum(Myotubes));
TotA=s1*s2;
MTDens=MTD/TotA;
end


function [xt1,yt1,xt2,yt2]=ManLab(Iin)
%title('Left click->add myoblast, Right click -> Remove myoblast, spacebar-> finish')
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
        figure(1);plot(xt1(j),yt1(j),'b.','MarkerSize',10)
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
    
    if button==3
        i = i+1;
        xt2(i) = xg; yt2(i)= yg;
        hold on
        figure(1);plot(xt2(i),yt2(i),'w.','MarkerSize',10)
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
end

end