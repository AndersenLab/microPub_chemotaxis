// lets you choose your folder
inputDir = getDirectory("Choose your input folder");

// the highest number you expect (95 in your case)
max = getNumber("How many images do you expect?", 6);

// lets you choose the naming of your files
SName = getString("Whats the strain prefix of your images?", "");

// makes a sub directory into your chosen folder
outputDir = inputDir + "\debug" + SName;
if(!File.exists(outputDir)) File.makeDirectory(outputDir);

// a loop through your images... (1 ... max)
for(i=1; i<=max; i++){
strain = inputDir + SName + i + ".jpg";

// opens the image
open(strain);

// draw quadrants (assay20180221)
setTool("line");
makeLine(1304, 43, 1304, 1903);
run("Add Selection...");
makeLine(380, 976, 2240, 975);
run("Add Selection...");
run("Set Scale...", "distance=1860 known=45.45 pixel=1 unit=mm global");
makeOval(1100, 770, 408, 408);
run("Add Selection...");

// Output measurements
run("Set Measurements...", "display add redirect=None decimal=3");

// setup multipoint tool
setTool("multipoint");
run("Point Tool...");

// Record measurments
run("Measure");

// Flatten image for saving
run("Flatten");

// saves the image into the sub directory as: overlay1..max
saveAs("jpg", outputDir + File.separator + "overlay" + i);

// closes all images before the next iteration
while(nImages > 0) close();
}

// dialogue box for user to export results
waitForUser("Export results", "Ensure results are correct for this strain.\nClick OK when done.");
selectWindow("Results");
saveAs("Measurements", inputDir + File.separator + SName + "Results.csv");

// Clear results for the next strain
run("Clear Results");
