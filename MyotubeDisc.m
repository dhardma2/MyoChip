filename= sprintf('stk_0026_pos7_FullMovie.czi kept stack.tif');
info=imfinfo(filename);
numframe=length(info);
%image size
flength=2203;
fwidth=2203;
f=1;
Raw(:,:,:,f) = imread(filename,f);
cooked=mat2gray(Raw(1:floor(fwidth),1:floor(flength),:,:));
%Find initial centroids
Iin=cooked(:,:,:,f);

Amax=[8000];
%Prefactors used to calibrate segmentation
x1=[1];
x2=[6];
x3=[2];
 [Myoblasts,Myotubes,NumMyoblasts,NumMyotubes,statsMB2,statsMT2]=Masks(Iin,x1,x2,x3,Amax,f);

%Select nuclei in myotubes

[xMT,yMT,Iin]=ManLab(Iin,Myotubes);
%draw ROIs
BWMBpatch=circledraw(flength,fwidth,xMT,yMT,Iin);
%remove myotube nuclei from image
Myotubes2=(Myotubes | BWMBpatch);

 
function [xout,yout,Iin]=ManLab(Iin,BWMB)
figure; hold on;
imshow(imoverlay(Iin,BWMB))
[sizex]=size(Iin(:,1,:,1));
[sizey]=size(Iin(1,:,:,1));
for row = 1 : 250 : sizex(1)%3154
 line([1, sizey(2)], [row, row]);
end
for column = 1 : 250 : sizey(2)%3954
 line([column , column ], [1, sizex(1)]);
end
x = 0;
y = 0;
button = 1;
i = 0;
while button <=3
    [xg,yg,button] = ginput(1);
    if button==1
        i = i+1;
        xout(i) = xg; yout(i)= yg;
        hold on
        plot(xout(i),yout(i),'b.','MarkerSize',25)
        drawnow
    end


end
end

function BWpatch=circledraw(length,width,xc,yc,Iin)
BWpatch=false(length,width);
for ii=1:max(size(xc))
    
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = length;
imageSizeY = width;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = xc(ii);
centerY = yc(ii);
radius = 15;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
%image(circlePixels) ;
colormap([0 0 0; 1 1 1]);
BWpatch= BWpatch | circlePixels;

end
end

function [Myoblasts,Myotubes,NumMyoblasts,NumMyotubes,statsMB2,statsMT2]=Masks(Iin2,x1,x2,x3,Amax,f)
BW4 = BinaryImage(Iin2,x1(1),x2(1),x3(1),2,2,Amax(1),300);
Myoblasts1=BW4(:,:,1);
Myotubes1=BW4(:,:,2);

se = strel('disk',20);
%pre-process using an 'Opening-by-reconstruction' method to define
%background/foreground.
Ie = imerode(Iin2,se);
Iobr = imreconstruct(Ie,Iin2);

BW5 = BinaryImage(Iobr,x1(1),x2(1),x3(1),4,2,Amax(1),300);
Myoblasts2=BW5(:,:,1);
Myotubes2=BW5(:,:,2);

%Combine pre-processed and original binary masks to provide optimum
%segmentation.
Myotubes=Myotubes1 | Myotubes2;
%figure; imshow(Myotubes)
Myoblasts=Myoblasts2 & ~Myotubes;

%figure; imshow(Myoblasts)

MaxXMB=size(Myoblasts);

%statsMT = regionprops('table',Myotubes,'Centroid','area','perimeter','boundingbox','MajorAxisLength');

%AreaMT=statsMT.Area;
%BB=statsMT.BoundingBox;

statsMB = regionprops('table',Myoblasts,'Centroid','area','perimeter','boundingbox','MajorAxisLength','Extrema'); 
%Function for moving large regions from myoblast to myotube mask
[Myoblasts, Myotubes] = movemyo(statsMB,Myoblasts,Myotubes);
Myoblasts=bwareaopen(Myoblasts,350);
statsMB2 = regionprops('table',Myoblasts,'Centroid','area','perimeter','boundingbox','MajorAxisLength');
statsMT2 = regionprops('table',Myotubes,'Centroid','area','perimeter','boundingbox','MajorAxisLength');
%statMB = regionprops(Myoblasts,Iin2,{'Centroid'});
%statMT = regionprops(Myotubes,Iin2,{'Centroid'});

numMB = numel(statsMB2.Area);
numMT = numel(statsMT2.Area);

NumMyoblasts(f)=numMB;
NumMyotubes(f)=numMT;

%figure;
%imshow(labeloverlay(Iin2,Myotubes))
%imshow(BW2)
end

function BW4 = BinaryImage(Iin,x1,x2,x3,strelmagd,strelmage,Amax,Amin)

            
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
    bigholes=bwareaopen(holes, 1000);
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
 
function [Myoblasts, Myotubes] = movemyo(statsMB,Myoblasts,Myotubes)
LengthMB=statsMB.MajorAxisLength;
%Move larger regions from Myoblasts to Myotubes mask
hold on
c=0;
for d=1:size(LengthMB)
    if LengthMB(d)>175
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
