clear
% Code to count Myocte and myotube nuclei 
%Import selected frames
prompt = {'Pixel width (microns):','Voxel depth:','Stack Size:','Mean Nucleus Volume:','Mean Myocyte Volume'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.2271617','0.5683904','100','2050','5000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%voxel size
vs=str2double(answer(1));
vd=str2double(answer(2));
voxvol=vs*vs*vd;
StackSize=str2double(answer(3));
MCNucMean = str2double(answer(4));
MCNucMean = round(MCNucMean / voxvol,0);
MCCellMean =str2double(answer(5));
MCCellMean = round(MCCellMean / voxvol,0);


for f=1:StackSize

RawA(:,:,:,f) = imread('Actin.tif',f);


RawN(:,:,:,f) = imread('Nuclei.tif',f);


end
cooked=mat2gray(RawA);
BWA=imbinarize(cooked(:,:,:));
seA = strel('diamond',1);
BWAdil = imdilate(BWA,seA);
BWfill = imfill(BWAdil,'holes');
BWAfilt=bwareaopen(BWfill,1000);
BWAMCfilt=bwareaopen(BWfill,MCCellMean);
CCA = bwconncomp(BWAMCfilt);
SA = regionprops(CCA,'Centroid','Area');
CentroidsA = cat(1,SA.Centroid);
AreaA = cat(1,SA.Area);
%figure;scatter3(CentroidsA(:,1),CentroidsA(:,2),CentroidsA(:,3))


cooked=mat2gray(RawN);
BW=imbinarize(cooked(:,:,:));
se = strel('diamond',1);
BWdil = imdilate(BW,se);
BWfill = imfill(BWdil,'holes');
BWfilt=bwareaopen(BWfill,round(0.75*MCNucMean,0));
CC = bwconncomp(BWfilt);
S = regionprops(CC,'Centroid','Area');
Centroids = cat(1,S.Centroid);
Area = cat(1,S.Area);
figure;scatter3(Centroids(:,1),Centroids(:,2),Centroids(:,3),'k')
hold on
%fN=isosurface(BWfilt);
%p=patch(fN,'FaceColor','blue','EdgeColor','none','FaceAlpha',.3);

%fA=isosurface(BWAMCfilt);
%pA=patch(fA,'FaceColor','red','EdgeColor','none','FaceAlpha',.3);

MCInd(1:length(Centroids))=0;
for k = 1 : length(Centroids)
    pointsImage = false(size(BW));
    Xval = round(Centroids(k,1));
    Yval = round(Centroids(k,2));
    Zval = round(Centroids(k,3));

    pointsImage(Yval, Xval, Zval) = true;
    MCCentroids = BWAMCfilt & pointsImage;
    if sum(sum(sum(MCCentroids))) ~= 0
        MCInd(k)=1;
    end
        
end


TNSing=0;
TNMult=0;
for i=1:length(Area)
    if MCInd(i) > 0 
        if Area(i)<2*MCNucMean
            TNSing=TNSing+1;
        else
            TNMult=TNMult+round((Area(i)/MCNucMean),0);
        end
    end
end
TotNuclei = TNSing + TNMult;