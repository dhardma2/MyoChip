for (i = 12; i < 22; i++) {
Expt="10 Hz";
Dir="D:/OptoImages_140922/";
setBatchMode(true);
if(i < 10){
Importer="open=["+Dir+Expt+"/"+Expt+"-01-Stitching-0"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
}
else{
Importer="open=["+Dir+Expt+"/"+Expt+"-01-Stitching-"+i+".czi] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2";
}	
run("Bio-Formats Importer",Importer );
rename("Pos"+i);
Pos=getTitle;
//Day=2;
run("Z Project...", "projection=[Average Intensity]");
run("Subtract Background...", "rolling=70 stack");
run("Split Channels");

selectWindow(getTitle);
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-MTDNA.tif");
setAutoThreshold("Default dark");
setThreshold(450, 66000);
setOption("BlackBackground", true);
run("Convert to Mask");
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-MTDNABW.tif");
close();

selectWindow(getTitle);
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-Actin.tif");
setAutoThreshold("Default dark");
setThreshold(1000, 66000);
setOption("BlackBackground", true);
run("Convert to Mask");
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-ActinBW.tif");
close();


selectWindow(getTitle);
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-BTX.tif");
close();

selectWindow(getTitle);
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-DNA.tif");
setAutoThreshold("Default dark");
setThreshold(1000, 66000);
//run("Threshold...");
setOption("BlackBackground", true);
run("Convert to Mask");              
saveAs("Tiff", Dir+Expt+"/"+Expt+"-"+Pos+"-DNABW.tif");             
close();

close();
};