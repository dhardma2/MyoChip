clear
prompt = {'Number of first video:','Number of last video:','Filename suffix:','Min branch length:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','3','s5_down_0','50'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
vid1=str2double(answer(1));
vid2=str2double(answer(2));
fname=char(answer(3));
MinBL=str2double(answer(4));

for v=1:(vid2-vid1+1)
    vidno=vid1+(v-1);
    filename= sprintf('%s%d.tif',fname,vidno);
    MTOriginal = imread(filename);
    MTBW = imbinarize(MTOriginal);
    Lmax=floor(length(MTBW)/3); 
    k=0;
    for i=1:3
        for j=1:3
            k=k+1;
            MTBWpatch=MTBW((i-1)*Lmax+1:i*Lmax,(j-1)*Lmax+1:j*Lmax);
            BWskelpatch = bwskel(MTBWpatch,'MinBranchLength',MinBL);
            BWBranchpatch = bwmorph(BWskelpatch,'branchpoints');
            BWBranchpatch = bwmorph(BWBranchpatch,'thicken',1);
            BWsplitpatch = BWskelpatch&~BWBranchpatch;

        %Remove very small regions (assumed artifacts)
            BWcleanpatch = bwareaopen(BWsplitpatch,MinBL);
            sPatch=regionprops(BWcleanpatch,'orientation','Area');
            AreaPatch = cat(1,sPatch.Area);
            AnglePatch=cat(1,sPatch.Orientation);
            Angrad = AnglePatch * pi /180;
            Angrad = 2 * Angrad; 
            AngleStats(k) = circ_stats(Angrad);
                        
        end
    end
    AngMeanDeg = cat(1,AngleStats.mean) *90 / pi;
    AngMedianDeg = cat(1,AngleStats.median) *90 / pi;
    AngSTDDeg=cat(1,AngleStats.std) *90 / pi;
    filename2=sprintf('%s%d_stats',fname,vidno);
    save(filename2,'AngleStats','AngMeanDeg','AngMedianDeg','AngSTDDeg');
    clear
end