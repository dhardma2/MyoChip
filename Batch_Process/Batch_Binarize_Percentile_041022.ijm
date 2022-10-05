for (i = 21; i < 22; i++) {
Expt="Trained";
Dir="E:/Optogen/";
Day="2";
//pixel size used for rolling ball background subtraction
SubtractionPxSize="70";
//apply batch mode to 
setBatchMode(true);
if(i < 10){
Importer="open=["+Dir+"Day"+Day+" "+Expt+"-01-Stitching-"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
}
else{
Importer="open=["+Dir+"Day"+Day+" "+Expt+"-01-Stitching-"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
}

//if(i < 10){
//Importer="open=["+Dir+Expt+"/Day"+Day+" Untrained "+Expt+"-01-Stitching-0"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
//}
//else{
//Importer="open=["+Dir+Expt+"/Day"+Day+" Untrained "+Expt+"-01-Stitching-"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
//}	
//import .czi images
run("Bio-Formats Importer",Importer );
rename("Pos"+i);
Pos=getTitle;
//z-stack based on average intensity
run("Z Project...", "projection=[Average Intensity]");
//background subtraction
run("Subtract Background...", "rolling="+SubtractionPxSize+" stack");
run("Split Channels");

function Percentile_Threshold(percentage){
	//percentage = 75; 
nBins = 256; 
resetMinAndMax(); 
getHistogram(values, counts, nBins); 
// find culmulative sum 
nPixels = 0; 
for (i = 0; i<counts.length; i++) 
  nPixels += counts[i]; 
nBelowThreshold = nPixels * percentage / 100; 
sum = 0; 
for (i = 0; i<counts.length; i++) { 
  sum = sum + counts[i]; 
  if (sum >= nBelowThreshold) { 
    setThreshold(values[0], values[i]); 
    print(values[0]+"-"+values[i]+": "+sum/nPixels*100+"%"); 
    i = 99999999;//break 
  } 
} 
}

function Make_Binary(im_title,percentage){
selectWindow(getTitle);
saveAs("Tiff", im_title+".tif");
Percentile_Threshold(percentage);
run("Convert to Mask");
run("Invert");
saveAs("Tiff", im_title+"BW.tif");             
close();

}

title=Dir+Expt+"/"+Expt+"-"+Pos+"-DNA";
Make_Binary(title,90);
title=Dir+Expt+"/"+Expt+"-"+Pos+"-BTX";
Make_Binary(title,75);
title=Dir+Expt+"/"+Expt+"-"+Pos+"-Actin";
Make_Binary(title,75);
title=Dir+Expt+"/"+Expt+"-"+Pos+"-MTDNA";
Make_Binary(title,95);

close();
};