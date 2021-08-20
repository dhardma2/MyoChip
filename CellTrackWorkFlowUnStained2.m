%Matlab code to segment, track and obtain metrics of motion and density for myocytes
clear;
nf=0;
%input video number for tracking and image parameters
prompt = {'Pixel width (microns):','Frame length (mins):','Number of first video:','Number of last video:','Filename suffix:','Position:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.65','5','2','6','_MTDay0-1 medium 1-1 5min','9'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%Track cells 'live' in GUI ('yes' increases duration of tracking)
answer2 = questdlg('Display live tracking?','Yes','No');

switch answer2
    case 'Yes'
        GUI = 1;
    case 'No'
        GUI = 0;
end

%Pixel size (microns)
ps=str2double(answer(1));
%Frame length (mins)
fl=str2double(answer(2));
vid1=str2double(answer(3));
vid2=str2double(answer(4));
suffix=char(answer(5));
position=str2double(answer(6));

for v=1:(vid2-vid1+1)
    if vid1 >= 10
        pref='00';
    else 
        pref='000';
    end
    vidno=vid1+(v-1);
filename= sprintf('stk_%s%d%s-pos%d.tif',pref,vidno,suffix,position);
info=imfinfo(filename);
numframe=6;%length(info);
%image size
flength=round(1000/ps,0);
fwidth=round(1000/ps,0);
for f=1:numframe
    nf=nf+1;
    g=(v-1)*numframe+f;
    %ff=f+90;
Raw(:,:,:,g) = imread(filename,f);


end
end
cooked=mat2gray(Raw(1:floor(fwidth),1:floor(flength),:,:));
%Find initial centroids
f=1;
Iin=cooked(:,:,:,f);

Amax=[8000];
%Prefactors used to calibrate segmentation
x1=[1];
x2=[6];
x3=[2];

[Myoblasts,Myotubes,NumMyoblasts,NumMyotubes,statsMB2,statsMT2]=Masks(Iin,x1,x2,x3,Amax,f,ps);

xt1(:,f)=statsMB2.Centroid(:,1);
yt1(:,f)=statsMB2.Centroid(:,2);
figure(1);imshow(Iin)
hold on
jj=length(xt1(:,1));
figure(1);plot(xt1(:,1), yt1(:,1), 'r.', 'MarkerSize', 10)
%a = [1:kk]'; b = num2str(a); c = cellstr(b);
%dx = 1; dy = 1; % displacement so the text does not overlay the data points
%text(xt1(:,1)+dx, yt1(:,1)+dy, c, 'Fontsize', 10);

%Function to manually label cells missed by segmentation 
[xt2,yt2,xb,yb]=ManLab(jj,xt1,yt1,Iin);
Myocytetot(f)=length(xt2)-length(xb);

kk=length(xt1(:,1));
Buff(1:kk)=0;

for f=2:nf
    f
Iin2=cooked(:,:,:,f);
[Myoblasts,Myotubes,NumMyoblasts,NumMyotubes,statsMB2]=Masks(Iin2,x1,x2,x3,Amax,f,ps);
numObj = numel(statsMB2.Area);

for ii=1:kk
    xt1(ii,f)=0;
    yt1(ii,f)=0;
    xmin=0; xmax=0;
    ymin=0; max=0;
    %initialise region of influence
    d=round(45.4/ps,0);
    g=f-1;
if xt1(ii,g) ~=0 && yt1(ii,g) ~=0  
    xmin=xt1(ii,g)-(d/2);
    xmax=xt1(ii,g)+d/2;
    ymin=yt1(ii,g)-(d/2);
    ymax=yt1(ii,g)+d/2;
        
     nc = 0;
     dist=10000;
     xt2=0;
     yt2=0;
     xnew=0;
     ynew=0;

for k=1:numObj
    xt2=statsMB2.Centroid(k,1);
    yt2=statsMB2.Centroid(k,2);
    
    
    if (xmin<xt2) && (xt2<xmax)
        
         if (ymin<yt2) && (yt2<ymax)
             nc=nc+1;
             xtnew(nc)=xt2;
             ytnew(nc)=yt2;
             index(nc)=k;
         end
    end
end
         
                if nc==0
                    xt1(ii,f)=0;  yt1(ii,f)=0;
                elseif nc>1
                for jj=1:nc
                    distnew=distanceA(xtnew(jj),xt1(ii,g),ytnew(jj),yt1(ii,g));   
                    if distnew<dist
                        xt2=xtnew(jj);
                        yt2=ytnew(jj);
                        dist=distnew;
                    end
                    xt1(ii,f)=xt2;
                    yt1(ii,f)=yt2;
                end
                else
                     xt1(ii,f)=xtnew(1);
                     yt1(ii,f)=ytnew(1);
                end 
              
    %Test for uniqueness
              for u=1:(ii-1)
                  if (xt1(ii,f)>0) && (yt1(ii,f)>0)
                    if (xt1(ii,f)==xt1(u,f)) &&  (yt1(ii,f)==yt1(u,f))
                     dist1=distanceA(xt1(ii,f),xt1(ii,g),yt1(ii,f),yt1(ii,g));
                     dist2=distanceA(xt1(u,f),xt1(u,g),yt1(u,f),yt1(u,g));
                         if dist1<dist2
                             xt1(u,f)=xt1(u,g); yt1(u,f)=yt1(u,g);
                             %figure(1);plot( xt1(u,g),  yt1(u,g), 'gx', 'MarkerSize', 10)
                         else
                             xt1(ii,f)=xt1(ii,g); yt1(ii,f)=yt1(ii,g);
                             %figure(1);plot( xt1(ii,g),  yt1(ii,g), 'gx', 'MarkerSize', 10)
                         end
                     
                    end
                  end
              end
end
end

    %plot cell motion 
    for ii=1:kk
       
     %Include a buffer in case segmentation has missed the
     %Myoblast in a given frame.
               if (xt1(ii,f)==0) &&  (yt1(ii,f)==0)
                  Buff(ii)=Buff(ii)+1;
                   if Buff(ii)< 3
                    xt1(ii,f)=xt1(ii,g);  yt1(ii,f)=yt1(ii,g);
                   end
               elseif xt1(ii,f)~=xt1(ii,g)
                   Buff(ii)=0;
                                    
               end
        if GUI == 1
            if  (xt1(ii,f)==0) && (yt1(ii,f)==0) && (xt1(ii,g)~=0)
                hold on
                figure(1);plot( xt1(ii,g),  yt1(ii,g), 'gx', 'MarkerSize', 10)
         
            else
            
                hold on
                figure(1);plot(xt1(ii,f),  yt1(ii,f), 'r.', 'MarkerSize', 20)
                plot([xt1(ii,g) xt1(ii,f)],[yt1(ii,g) yt1(ii,f)],'r');
          
            end
       
        end
     end  
 
                
end
%*************************************************************************************************
%Extracting metrics from tracking data

%Extract cells which have not been removed prematurely (more than 5
%recorded Frames).
%Removes the effects of the buffer used to placehold cells if they are 
%missed by segmentation in a frame.
[celldatx,celldaty,celldatind,D0]=ExtDat(xt1,yt1);

%Category 1: Discrete averages

%Calculate the density of myoblasts /mm^2

s1=flength;
s2=fwidth;
Densum=length(celldatx)/(s1*ps*s2*ps);
Densmm=Densum*1000*1000;

%Calculate the average nearest neighbour distance between cells and the 
%'nearest neighbour value' (Clustered: NNval->0, Randomly distributed->1,
%uniform distribution->2.15)
[NNd, SNNd, Zval, ZvalBC]=MeanNearestNeighbour(s1,s2,celldatx,celldaty,ps);
MeanNNd=NNd;
SDNNd=SNNd;
ZvalNNd=Zval;
ZvalNNdBC=ZvalBC;

%Calculating the normalised average Voronoi Cell area
[MeanVAnorm, SDVAnorm]=VorArea(celldatx,celldaty,s1,s2,ps);
MeanVA=MeanVAnorm;
SDVA=SDVAnorm;

%Category 2

%Calculate total distance travelled by each myoblast and find the mean
%speed.
[totd, totsd, meanspeed]=celvel(celldatx,celldaty,celldatind);
%Convert units into microns and microns per minute
totdc=totd*ps;
meanspeedc=meanspeed*(ps/fl);

%Calculate the resultant vectors, persistence, velocity and angle of
%motion.
[ResVectD,Persistence, meanvel,theta,weight]=VectAng(celldatx,celldaty,celldatind,totd);
%Convert units into microns and microns per minute
ResVectDc=ResVectD*(ps);
meanvelc=meanvel*(ps/fl);

%Calculate the magnitude, displacement and persistence of angular motion.
[rotsp,rotvel,rotpers]=AngVel(celldatx,celldaty,celldatind,fl);
t=1;
DensityMBFinal(t)=mean(Densmm);
NearestNeighbourFinal(t)=mean(MeanNNd);
NearestNeighbourSDFinal(t)=mean(SDNNd);
NearestNeighbourZvalFinal(t)=mean(ZvalNNdBC);
NearestNeighbourHyp(t)=0;
if NearestNeighbourZvalFinal(t)<=-1.96
    NearestNeighbourHyp(t)=-1;
elseif NearestNeighbourZvalFinal(t)>=1.96
    NearestNeighbourHyp(t)=1;
end    
VoronoiAreaFinal(t)=mean(MeanVA);
VoronoiAreaSDFinal(t)=mean(SDVA);
MeanMBSpeed(t)=mean(meanspeedc);
SDMBSpeed(t)=std(meanspeedc);
MeanMBVel(t)=mean(meanvelc);
SDMBVel(t)=std(meanvelc);
MeanMBPersistence(t)=mean(Persistence);
SDMBPersistence(t)=std(Persistence);
MeanMBRotSpeed(t)=mean(rotsp);
SDMBRotSpeed(t)=std(rotsp);
MeanMBRotVel(t)=mean(rotvel);
SDMBRotVel(t)=std(rotvel);
MeanMBRotPers(t)=mean(rotpers);
SDMBRotPers(t)=std(rotpers);

%Run Rayleigh-Moore test for Vector Angles
[r_s, resultant_phase] = moore_raleigh(theta, ResVectDc);
r_star(t)=r_s;
if r_star(t)>0.999
    RMTestMB(t)=0;
else
    RMTestMB(t)=1;
end
%********************************************************************************
%clear large matrices for saving
clear cooked Iin Iin2

   
%Create binary mask for Myoblasts and Myotubes
function [Myoblasts,Myotubes,NumMyoblasts,NumMyotubes,statsMB2,statsMT2]=Masks(Iin2,x1,x2,x3,Amax,f,ps)
%critical length of myocyte before it is classed as a myotube
crit_myoc_len=round(79.45/ps,0);
Amin=round(136.2/ps,0);
BW4 = BinaryImage(Iin2,x1(1),x2(1),x3(1),2,2,Amax(1),Amin,ps);
Myoblasts1=BW4(:,:,1);
Myotubes1=BW4(:,:,2);
%create structuring element (disk shaped)
se = strel('disk',round(9.08/ps));
%pre-process using an 'Opening-by-reconstruction' method to define
%background/foreground.
Ie = imerode(Iin2,se);
Iobr = imreconstruct(Ie,Iin2);

BW5 = BinaryImage(Iobr,x1(1),x2(1),x3(1),4,2,Amax(1),Amin,ps);
Myoblasts2=BW5(:,:,1);
Myotubes2=BW5(:,:,2);

%Combine pre-processed and original binary masks to provide optimum
%segmentation.
Myotubes=Myotubes1 | Myotubes2;

Myoblasts=Myoblasts2 & ~Myotubes;

MaxXMB=size(Myoblasts);

statsMB = regionprops('table',Myoblasts,'Centroid','area','perimeter','boundingbox','MajorAxisLength','Extrema'); 
%Function for moving large regions from myoblast to myotube mask
[Myoblasts, Myotubes] = movemyo(statsMB,Myoblasts,Myotubes,crit_myoc_len);
Myoblasts=bwareaopen(Myoblasts,round(158.9/ps,0));
statsMB2 = regionprops('table',Myoblasts,'Centroid','area','perimeter','boundingbox','MajorAxisLength');
statsMT2 = regionprops('table',Myotubes,'Centroid','area','perimeter','boundingbox','MajorAxisLength');

numMB = numel(statsMB2.Area);
numMT = numel(statsMT2.Area);

NumMyoblasts(f)=numMB;
NumMyotubes(f)=numMT;

end

function BW4 = BinaryImage(Iin,x1,x2,x3,strelmagd,strelmage,Amax,Amin,ps)

            
I=Iin(:,:);

%Changes in contrast can be detected by operators that calculate the gradient of an image. 
%To create a binary mask containing the segmented cell, calculate the gradient image and apply a threshold.

%Use edge and the Sobel operator to calculate the threshold value. 
%Tune the threshold value and use edge again to obtain a binary mask that contains the segmented cell.

[~,threshold] = edge(I,'sobel');

fudgeFactor = 0.5;
BWs = edge(I,'sobel',threshold * fudgeFactor);
%Display the resulting binary gradient mask.
%imshow(BWs)

%The binary gradient mask shows lines of high contrast in the image. 
%These lines do not quite delineate the outline of the object of interest. 
%Compared to the original image, there are gaps in the lines surrounding the object in the gradient mask. 
%These linear gaps will disappear if the Sobel image is dilated using linear structuring elements.
%Create two perpindicular linear structuring elements by using strel function.
sm=strelmagd;
se90 = strel('line',sm,90);
se0 = strel('line',sm,0);
%Dilate the binary gradient mask using the vertical structuring element followed by the horizontal structuring element. 
%The imdilate function dilates the image.
BWsdil = imdilate(BWs,[se90 se0]);
%imshow(BWsdil)
%The dilated gradient mask shows the outline of the cell quite nicely, but there are still holes in the interior of the cell. 
%To fill these holes, use the imfill function.

%The dilated gradient mask shows the outline of the cell quite nicely, but there are still holes in the interior of the cell. 
%To fill these holes, use the imfill function.
    
    BWdfill = imfill(BWsdil,'holes');
    %Myotubes can come together to make large holes. These must be
    %discounted from the fill.
    %figure; imshow(BWdfill)
    holes=BWdfill & ~BWsdil;
    %figure; imshow(holes)
    bigholes=bwareaopen(holes, round(454/ps,0));
    %figure; imshow(bigholes)
    smallholes = holes & ~bigholes;
    %figure; imshow(smallholes)
    BWdfill = BWsdil | smallholes;
    %figure; imshow(labeloverlay(Iin,BWdfill))
%imshow(BWdfill)

sd2=strelmage;
seD = strel('diamond',sd2);
BW2 = imerode(BWdfill,seD);
BW2 = imerode(BW2,seD);
%imshow(BW2)

%remove regions greater than Amax pixels in area

BWAgg2=bwareaopen(BW2,Amax);


seA=strel('disk',x1*sd2);

BWAgg2a=imerode(BWAgg2,seA);

BWAgg2=BWAgg2a-bwareaopen(BWAgg2a,Amax);


Myotubes=BWAgg2a-BWAgg2;


%Filter larger regions
BW2=BW2-bwareaopen(BW2,Amax);
BW2=bwareaopen(BW2,Amin);

sizex=size(BWAgg2(:,1));
sizey=size(BWAgg2(1,:));
for i=1:sizex
    for j=1:sizey(2)
        if (BW2(i,j)==0) && (BWAgg2(i,j)==1)
            BW2(i,j)=1;
        end
    end
end


%Filter smaller regions
BWAgg1=BW2-bwareaopen(BW2,x2*Amin);

seA=strel('disk',x3*sm);

BWAgg1=imdilate(BWAgg1,seA);
%imshow(BWAgg1)

for i=1:sizex
    for j=1:sizey(2)
        if (BW2(i,j)==0) && (BWAgg1(i,j)==1)
            BW2(i,j)=1;
        end
    end
end


%remove areas less than 10 pixels in area
BW2=bwareaopen(BW2,Amin);
%figure
%imshow(BW2)
BW4(:,:,1)=BW2;
BW4(:,:,2)=Myotubes;

end
 
function [Myoblasts, Myotubes] = movemyo(statsMB,Myoblasts,Myotubes,crit_myoc_len)
LengthMB=statsMB.MajorAxisLength;
%Move larger regions from Myoblasts to Myotubes mask
hold on
c=0;
for d=1:size(LengthMB)
    if LengthMB(d)>crit_myoc_len
        c=c+1;
       
        Extrema=cell2mat(statsMB.Extrema(d));
        %plot(xC,yC,'r.','MarkerSize',20)
        x=ceil(Extrema(1,1));
        y=ceil(Extrema(1,2));
        i=0;
        while Myoblasts(y,x)==false
            x=x+1;
            y=y+1;
            i=i+1;
            if x>MaxXMB(2)
                
                x=MaxXMB(2);
            end
            if y>MaxXMB(1)
                
                y=MaxXMB(1);
            end
            if i>100
                break
            end
        end
           % plot(x,y,'r.','MarkerSize',10)
               if Myoblasts(y,x)==true
                 marker = false(size(Myoblasts));
                 marker(y,x) = true;
                 im = imreconstruct(marker,Myoblasts);
                 %figure;
                 %imshow(im)
                 Myoblasts=Myoblasts & ~im;
                 Myotubes=Myotubes | im;
               end
       
        
     end
end

end

function [xt1,yt1,xb,yb]=ManLab(kk,xt1,yt1,Iin)
[sizex]=size(Iin(:,1,:,1));
[sizey]=size(Iin(1,:,:,1));
for row = 1 : 500 : sizex(1)%3154
 line([1, sizey(2)], [row, row]);
end
for column = 1 : 500 : sizey(2)%3954
 line([column , column ], [1, sizex(1)]);
end
x = 0;
y = 0;
button = 1;
i = 0;
j = 0;
while button <=3
    [xg,yg,button] = ginput(1);

    if button==1
        j = j+1;
        xt1(j,1) = xg; yt1(j,1)= yg;
        hold on
        figure(1);plot(xt1(j),yt1(j),'r.','MarkerSize',10)
        %a = [j]'; b = num2str(a); c = cellstr(b);
        %dx = 1; dy = 1; % displacement so the text does not overlay the data points
        %text(xt1(j,1)+dx, yt1(j,1)+dy, c, 'Fontsize', 10);
        drawnow
        
    end
    if button==3
      i = i+1;
        xb(i,1) = xg; yb(i,1)= yg;
        hold on
        plot(xb(i),yb(i),'g.','MarkerSize',20)
        drawnow
    end
end
end

function distnew=distanceA(x2,x1,y2,y1)
    distnew=sqrt(((x2-x1)^2)+((y2-y1)^2));
end

function [celldatx,celldaty,celldatind,D0]=ExtDat(xt1,yt1)
ind=0;
D0=0;
for i=1:size(xt1,1)
    t=0;
    
 for j=1:size(xt1,2)
     if xt1(i,j)~=0 && yt1(i,j)~=0
         t=t+1;
     end
 end
 if xt1(i,t)==xt1(i,1) && yt1(i,t)==yt1(i,1)
     D0=D0+1;
 elseif t>5 
     ind=ind+1;
     celldatind(ind)=i;
     celldatx(ind,:)=xt1(i,:);
     celldaty(ind,:)=yt1(i,:);
 end
end

%Remove buffer effect
s=size(celldatx,2);
for k=1:ind
    for m=1:s
        if celldatx(k,s+1-m)==0 && celldatx(k,s-m) ~=0
            celldatx(k,s-m)=0;
            celldatx(k,s-1-m)=0;
            celldaty(k,s-m)=0;
            celldaty(k,s-1-m)=0;
            break
        end
            
    end
end

end

function [MeanNNd,SDNNd,ZvalNNd,ZvalNNdBC]=MeanNearestNeighbour(s1,s2,celldatx,celldaty,ps)
xy(1,:)=celldatx(:,1);
xy(2,:)=celldaty(:,1);
[idx] = nearestneighbour(xy);
L=length(idx);
for i=1:L
    x1=xy(1,i);
    y1=xy(2,i);
    x2=xy(1,idx(i));
    y2=xy(2,idx(i));
    dx=x2-x1;
    dy=y2-y1;
    NNdist(i)=sqrt(dx*dx+dy*dy);    
end
MeanNNd=sum(NNdist(:))/L;
SampleArea=s1*s2;
SamplePerim=2*(s1+s2);
RefNNDist=0.5*sqrt(SampleArea/L);
%Boundary effect correction
BC=(0.0514+(0.041/sqrt(L)))*(SamplePerim/L);
ExpNNd=RefNNDist+BC;
VarNNd=0.07*(SampleArea/(L*L));
VarNNdBC=0.07*(SampleArea/(L*L))+(0.037*SamplePerim*sqrt(SampleArea/L^5));
ZvalNNd=(MeanNNd-ExpNNd)/sqrt(VarNNd);
ZvalNNdBC=(MeanNNd-ExpNNd)/sqrt(VarNNdBC);

MeanNNd=MeanNNd*ps;
SDNNd=std(NNdist)*ps;
end

function [MeanVAnorm, SDVAnorm]=VorArea(celldatx,celldaty,s1,s2,ps)
xya(:,1)=celldatx(:,1);
xya(:,2)=celldaty(:,1);
xy=unique(xya,'rows');
figure; [v,c]=voronoin(xy);
xlim([0 s2]); ylim([0 s1]);

w=0;
for i=1:length(c)
    cells=c{i};
    k=0;
    for j=1:length(cells)
        if v(cells(j),1)>0 && v(cells(j),2)>0 && v(cells(j),1)<s2 && v(cells(j),2)<s1
            k=k+0;
        else
            k=k+1;
        end
    end
        if k==0
            w=w+1;
            c2(w)=i;
            vorcell(w)=cell(c(i));
            
        end
end
vorcell=vorcell';
 
A = zeros(length(vorcell),1) ;
for i = 1:length(vorcell)
    v1 = v(vorcell{i},1) ; 
    v2 = v(vorcell{i},2) ;
    patch(v1,v2,rand(1,3))
    A(i) = polyarea(v1,v2) ;
    A(i)=A(i)*ps*ps;
end
MeanVA=mean(A);
MeanVAnorm=MeanVA/(s1*s2*ps*ps);
SDVA=std(A);
SDVAnorm=SDVA/(s1*s2*ps*ps);
end

function [totd, totsd, meanvel]=celvel(celldatx,celldaty,ind)
s=size(celldatx,2);
for k=1:size(ind,2)
[td,sd,mv]=distanceB(s,celldatx(k,:),celldaty(k,:));
totsd(k)=sd;
totd(k)=td;
meanvel(k)=mv;
end
end

function [ResVectD,Persistence,mv,theta,weight]=VectAng(celldatx,celldaty,ind,totd)
%Find dx
%find dy
%Use functions to find Resultant Vector and Theta
s=size(celldatx,2);
for k=1:size(ind,2)
    x1=celldatx(k,1); y1=celldaty(k,1);
    if celldatx(k,s)~=0
      x2=celldatx(k,s); y2=celldaty(k,s);
      weight(k)=s;
    else
        for m=2:s
        if celldatx(k,m)==0
            x2=celldatx(k,m-1);
            y2=celldaty(k,m-1);
            weight(k)=(m-1);
            break
        end
        end
    end
  
    dx=x2-x1;
    dy=y2-y1;
    [rvd]=ResVect(dx,dy);
    ResVectD(k)=rvd;
    Persistence(k)=ResVectD(k)/totd(k);
    mv(k)=ResVectD(k)/weight(k);
    [thet]=Ang(dx,dy);
    theta(k)=thet;
end
end

function [rotsp,rotvel,rotpers]=AngVel(celldatx,celldaty,ind,fl)
%Find dx
%find dy
%Use functions to find the magnitude and displacement of angular motion.
s=size(celldatx,2);
for k=1:size(ind,2)
    x1=celldatx(k,1); y1=celldaty(k,1);
    x2=celldatx(k,2); y2=celldaty(k,2);
    n=0;
    dx=x2-x1;
    dy=y2-y1;
    
    [thet]=Ang(dx,dy);
    theta1=thet;
    for kk=3:s
        
        if celldatx(k,kk)~=0
            n=n+1;
            x3=celldatx(k,kk); y3=celldaty(k,kk);
            dx=x3-x2;
            dy=y3-y2;
            [thet]=Ang(dx,dy);
            dtheta(n)=angdiff(theta1,thet);
            theta1=thet;
            x2=x3; y2=y3;
        end
    end
    dthetad=dtheta*(180/pi);
    rotvel(k)=sum(dthetad)/n;
    rotvel(k)=rotvel(k)/fl;
    rotsp(k)=0;
    for mm= 1:n
        rotsp(k)=rotsp(k)+sqrt(dthetad(mm)^2);
    end
    rotsp(k)=rotsp(k)/n;
    rotsp(k)=rotsp(k)/fl;
    if rotsp(k)==0
        rotspeed(k);
    end
    rotpers(k)=sqrt(rotvel(k)^2)/rotsp(k);
    
end
      
end

function [td,sd,ms]=distanceB(s,celldx,celldy)
n=0;
    for m=2:s
        if celldx(m)~=0
            n=n+1;
            d(m-1)=sqrt((celldx(m)-celldx(m-1))^2+(celldy(m)-celldy(m-1))^2);
        end
    end
    td=sum(d);
    sd=std(d);
    ms=td/n;

end

function [rvd]=ResVect(dx, dy)
rvd=sqrt((dx^2)+(dy^2));
end

function [thet]=Ang(dx,dy)
if (dx>=0) && (dy==0) 
    thet=0;
    
elseif (dx>0) && (dy>0) 
    thet=atan(dy/dx);
    
elseif (dx==0) && (dy>0) 
    thet=pi/2;

elseif (dx<0) && (dy>0) 
    thet=atan(dy/dx)+pi;

elseif (dx<0) && (dy==0) 
    thet=pi;
    
elseif (dx<0) && (dy<0) 
    thet=atan(dy/dx)+pi;
    
elseif (dx==0) && (dy<0)
    thet=1.5*pi;
    
elseif (dx>0) && (dy<0)
     thet=atan(dy/dx)+(2*pi);
end
end
