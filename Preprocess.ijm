rename("Pos3");
Pos=getTitle;
Day="D8"
dir = getDirectory("Choose a Directory ");
run("Z Project...", "projection=[Average Intensity]");
run("Subtract Background...", "rolling=70 stack");
run("Split Channels");

selectWindow(getTitle);
saveAs("Tiff", dir+Pos+"-HISPHD-"+Day+"-DNA.tif");              
close();
selectWindow(getTitle);
saveAs("Tiff", dir+Pos+"-HISPHD-"+Day+"-Phall.tif");
close();
selectWindow(getTitle);
saveAs("Tiff", dir+Pos+"-HISPHD-"+Day+"-MTDNA.tif");
close();
selectWindow(getTitle);
saveAs("Tiff", dir+Pos+"-HISPHD-"+Day+"-Actin.tif");
close();