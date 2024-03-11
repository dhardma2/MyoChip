%% Matlab code for discretising overlapping myotubes by angle of orientation
%% and providing a cell count. David Hardman, University of Edinburgh, 2024.
[file,path]=uigetfile;
filename=fullfile(path,file);
prompt = {'Min. branch length:','Min. cell length:','Bin angle:','Max image hole:','Binary dilation factor:','Display images(1/0):'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'50','500','30','2000','5','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%Minimum length of cell section to use when binarising image.
Min_branch_length = str2double(answer(1));
%Minimum cell length. Used to remove orphaned segments in post-processing.
Min_cell_length= str2double(answer(2));
%Bin angle. Smaller angle gives a greater chance of discretising cells
%which overlap with an acute angle but will mean less chance of recombining
%a cell which is curved.
D_Angle=str2double(answer(3));
N_sections=180/D_Angle;
%Maximum image hole size for preprocessing binary image
Max_hole=str2double(answer(4));
%Dilation factor for removing fuzz etc at image border
DilF=str2double(answer(5));

DispIm=str2double(answer(6));

image=imread(filename);
Im_size=size(image);
len=Im_size(1);
wid=Im_size(2);
BWAngsort=zeros(len,wid,N_sections);
%Create binary image.
MTBW=imbinarize(image);
%% Preprocessing: Cleaning the image for skeletonising effectively. 
%This section can be edited or removed depending on the original image quality.

%Remove noise and small sections.
BWCleaned=bwareaopen(MTBW,Max_hole);
%Fill holes
BWCleaned=~bwareaopen(~BWCleaned, Max_hole);
%Dilate to remove fuzzt borders.
seMT = strel('disk',DilF);
BWCleaned=imdilate(BWCleaned,seMT);

%% Skeletonise image and split the branch points
BWSkel=bwskel(BWCleaned,'MinBranchLength',Min_branch_length);
BWBranch=bwmorph(BWSkel,'branchpoints');
BWBranch=bwmorph(BWBranch, 'thicken',1);
BWSplit=BWSkel&~BWBranch;
Region_properties=regionprops(BWSplit,'orientation','PixelList');

%% Discretise cell segments by orientation angle
for i=1:length(Region_properties)
    Pixel_region=Region_properties(i).PixelList;
    for j=1:N_sections
        if Region_properties(i).Orientation+90 >(j-1)*D_Angle && Region_properties(i).Orientation+90 < j*D_Angle 
            for ii=1:length(Pixel_region)
                BWAngsort(Pixel_region(ii,2),Pixel_region(ii,1),j)=1;
            end
        end
    end
end

%% Dilate sections to re-join cells
BWAngsort=imdilate(BWAngsort,seMT);

%% Count cells
for k=1:N_sections
    BWAng=imbinarize(BWAngsort(:,:,k));
    BWAng=bwskel(BWAng);
    BWAng=bwareaopen(BWAng,Min_cell_length);
    BWAng=imdilate(BWAng,seMT);
    Angle_region_props=regionprops(BWAng);
    if DispIm==1
        figure;imshow(BWAng);
    end
    MT_count(k)=length(Angle_region_props);
end

Total_cells=sum(MT_count);
         
   
    
    