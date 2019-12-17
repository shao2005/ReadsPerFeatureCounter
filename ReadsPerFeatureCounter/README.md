# ReadsPerFeatureCounter

Requirements:

    •	Picard 2.6.0
    •	Text file with genes transcription in the format:
          <gene name>	<chromozom>	<start of exons in the transcript>	<end of exons in the transcript>
          Each line represents a different transcript of the gene. Therefor each gene can have more than 1 line.
          For example:
          g33123	chr17523	2529197,2533009,2539730	2529300,2533045,2539780
			
The program has 2 different mode: single cell and population.

Command for single cell run:

    java -jar ReadsPerFeatureCounter.jar sc  <pathToWorkingDirectory>  <BamFileName>  <FolderOfTheBamFile>  <pathToTranscriptsFile>  <outputFolder>  <geneToStart>  <uniqueId> <CellStart>  <cellEnd>
 
Command for population run:

    java -jar ReadsPerFeatureCounter.jar pop  <pathToWorkingDirectory>  <pathToBamFiles> <pathToTranscriptsFile>  <outputFolder>  <geneToStart>  <uniquId>  <FileName1>  <FileName2> <FileName3>…

cellStart/cellEnd – represent the cell interval we want to check. In case the interval is not continuous, the program can get more than one interval. For example, if the cells interval is 10-70 and 82-90, we send as parameters "10 70 82 90".

geneToStart – in case we want to start our process from specific gene. Otherwise '-'

unideID – 'u' if we want to count every pair of reads with the same id once. 'n' otherwise.

The program output separated into 6 files, written to outputFolder:

    1)	"unknownJunctions.csv" – calls that fall on unknown junctions.
    
    2)	"junctions.csv" – calls that fall on known junctions.
    
    3)	"introns.csv" – calls that fall on introns.
    
    4)	"exons.csv" – calls that fall on exons.
    
    5)	"geneCalls.csv" – all the calls for each gene.
    
    6)	"difCigar.csv" – calls with the cigar "*".

