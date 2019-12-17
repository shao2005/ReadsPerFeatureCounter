import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

import javax.xml.transform.sax.SAXTransformerFactory;

public class Gene{
	
	private String name;
	private String chrom;	
	protected GenePart[] exons;
	protected GenePart[] junctions;
	protected GenePart[] introns;
	private LinkedList<String> exonStarts;
	private LinkedList<String> exonEnds;
	private int[][] exonsArray;
	private int[][] intronsArray;	
	private List<int[]> unknowncells;
	private List<int[]> cells;
	private List<int[]> difCigarCells;
	private int[] unknownBF;
	private int[] totalCallsBF;
	private int[] difCigarBF;

	//constructors
	
	public Gene(String name, String chrom, GenePart[] exons, GenePart[] junctions,GenePart[] introns, LinkedList<String> exonStarts, LinkedList<String> exonEnds,
			List<int[]> unknowncells, List<int[]> cells, List<int[]> difCigarCells, int numOfBamFiles) {
		this.name = name;
		this.chrom = chrom;
		this.exons = exons;
		this.junctions = junctions;
		this.introns = introns;
		this.exonStarts = exonStarts;
		this.exonEnds = exonEnds;
		this.unknowncells = unknowncells;
		this.cells = cells;
		this.difCigarCells = difCigarCells;
		this.unknownBF=new int[numOfBamFiles];
		this.totalCallsBF=new int[numOfBamFiles];
		this.difCigarBF=new int[numOfBamFiles];

	}
		
	public Gene(String name, String chrom, List<cellsGroup> cells , LinkedList<String> exonStarts, LinkedList<String> exonEnds, int numOfBamFiles){
		this.name=name;
		this.chrom=chrom;
		this.exonStarts=exonStarts;
		this.exonEnds=exonEnds;
		if(cells!=null){
			this.unknowncells=cellsGroup.createCellsList(cells);
			this.cells=cellsGroup.createCellsList(cells);
			this.difCigarCells=cellsGroup.createCellsList(cells);
		}
		this.unknownBF=new int[numOfBamFiles];
		this.totalCallsBF=new int[numOfBamFiles];
		this.difCigarBF=new int[numOfBamFiles];
	}

//************************getters and setters***********************//
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public GenePart[] getExons() {
		return exons;
	}

	public void setExons(GenePart[] exons) {
		this.exons = exons;
	}

	public GenePart[] getJunctions() {
		return junctions;
	}

	public void setJunctions(GenePart[] junctions) {
		this.junctions = junctions;
	}

	public GenePart[] getIntrons() {
		return introns;
	}

	public void setIntrons(GenePart[] introns) {
		this.introns = introns;
	}

	public LinkedList<String> getExonStarts() {
		return exonStarts;
	}

	public void setExonStarts(LinkedList<String> exonStarts) {
		this.exonStarts = exonStarts;
	}

	public LinkedList<String> getExonEnds() {
		return exonEnds;
	}

	public void setExonEnds(LinkedList<String> exonEnds) {
		this.exonEnds = exonEnds;
	}

	public int[][] getExonsArray() {
		return exonsArray;
	}

	public void setExonsArray(int[][] exonsArray) {
		this.exonsArray = exonsArray;
	}

	public int[][] getIntronsArray() {
		return intronsArray;
	}

	public void setIntronsArray(int[][] intronsArray) {
		this.intronsArray = intronsArray;
	}

	
	public List<int[]> getCells(){
		return cells;
	}
	
	public void setCells(List<int[]> cells){
		this.cells=cells;
	}
	
	public List<int[]> getUnknowncells() {
		return unknowncells;
	}

	public void setUnknowncells(List<int[]> unknowncells) {
		this.unknowncells = unknowncells;
	}
	public List<int[]> getDifCigarCells() {
		return difCigarCells;
	}

	public void setDifCigarCells(List<int[]> difCigarCells) {
		this.difCigarCells = difCigarCells;
	}

	public int[] getUnknownBF() {
		return unknownBF;
	}

	public void setUnknownBF(int[] unknownBF) {
		this.unknownBF = unknownBF;
	}

	public int[] getTotalCallsBF() {
		return totalCallsBF;
	}

	public void setTotalCallsBF(int[] totalCallsBF) {
		this.totalCallsBF = totalCallsBF;
	}

	public int[] getDifCigarBF() {
		return difCigarBF;
	}

	public void setDifCigarBF(int[] difCigarBF) {
		this.difCigarBF = difCigarBF;
	}
	
	//*******************getters and setters end********************//

	
	//set exons array according to the starts and ends list
		public void setExons(List<List<String>> starts, List<List<String>> ends,List<cellsGroup> cells, int numOfBamFiles, char flag) {
			List<String> exonsString=new ArrayList<>();
			List<GenePart> exonsList=new ArrayList<>();
			//create a list of the exons
			for(int i = 0; i < starts.size(); i++){
				for(int j = 0; j < starts.get(i).size(); j++){			
					String exon = starts.get(i).get(j) + "-" + ends.get(i).get(j);
					if(!exonsString.contains(exon)){
						Exon e;
						if(flag == 's' || flag =='l')
							e = new Exon("", Util.convertToInteger(starts.get(i).get(j)), Util.convertToInteger(ends.get(i).get(j)), cellsGroup.createCellsList(cells), i, numOfBamFiles);
						else
							e = new Exon("", Util.convertToInteger(starts.get(i).get(j)), Util.convertToInteger(ends.get(i).get(j)), null, i, numOfBamFiles);

						exonsString.add(exon);
						exonsList.add(e);
					}
				}
			}
			//sort list
			Collections.sort(exonsList, new GeneComperator());
		
			//copy the list to an array
			exons = new Exon[exonsList.size()];
			for(int i = 0; i < exonsList.size(); i++)
				exons[i] = exonsList.get(i);
				
		}
		
		//set the junction array according to the introns array
		public void setJunctions(List<cellsGroup> cells, int numOfBamFiles, char flag) {
			junctions=new GenePart[introns.length]; 
			for(int i = 0; i < introns.length; i++){
				if(flag == 's' || flag == 'l')
					junctions[i] = new Junction("", introns[i].getStart(), introns[i].getEnd(),cellsGroup.createCellsList(cells), introns[i].getEnd() - introns[i].getStart(), i, 0, numOfBamFiles);//TODO: change the indexes
				else
					junctions[i] = new Junction("", introns[i].getStart(), introns[i].getEnd(),null, introns[i].getEnd() - introns[i].getStart(), i, 0, numOfBamFiles);//TODO: change the indexes

			}	
		}

		//set introns array according to the starts and ends list
		public void setIntrons(List<List<String>> starts, List<List<String>> ends,List<cellsGroup> cells, int numOfBamFiles, char flag) {
			List<GenePart> intronsList = new LinkedList<GenePart>();
			List<String> intronsString = new ArrayList<String>();
			
			//create a list of the introns
			for(int i = 0; i < starts.size(); i++){
				for(int j = 0; j < starts.get(i).size() - 1; j++){
					String intronName = ends.get(i).get(j) + "-" + starts.get(i).get(j+1);
					if(!intronsString.contains(intronName)){
						Intron intron;
						if(flag == 's' || flag =='l')
							intron = new Intron("", Util.convertToInteger(ends.get(i).get(j))+1, Util.convertToInteger(starts.get(i).get(j+1))-1, cellsGroup.createCellsList(cells), numOfBamFiles);
						else
							intron = new Intron("", Util.convertToInteger(ends.get(i).get(j))+1, Util.convertToInteger(starts.get(i).get(j+1))-1, null, numOfBamFiles);

						intronsString.add(intronName);
						intronsList.add(intron);
					}
				}
			}
			
			//sort list
			Collections.sort(intronsList, new GeneComperator());
			
			//copy the list to an array
			introns=new Intron[intronsList.size()];
			for(int i = 0; i < intronsList.size(); i++)
				introns[i] = intronsList.get(i);

		}
			

		//finds the start index of the junction
		public int findIndex(String cigar){
			int index = -1;
			if(cigar.contains("N"))
				index = cigar.indexOf((int)'N');
			if(cigar.contains("D") && (cigar.indexOf((int)'D')<index || index == -1) )
				index = cigar.indexOf((int)'D');
			if(cigar.contains("I")&& (cigar.indexOf((int)'I')<index || index == -1) )
				index = cigar.indexOf((int)'I');
			if(cigar.contains("S")&& (cigar.indexOf((int)'S')<index || index == -1) )
				index = cigar.indexOf((int)'S');
			return index;
		}
		
		//find the maximal size of the gene transcripts
		public int maxTrenscriptLength(List<List<String>> list){
			int maxLength = 0;
			for(int i = 0; i < list.size(); i++)
			{
				if(maxLength < list.get(i).size())
					maxLength = list.get(i).size();
			}
			return maxLength;
		}
		




		//add 1 junctionCall to specific cell
		//cellRange and cell - for singleCell, bamFieIndex - for population
		public boolean addToJunction(int start, int end, String ID, int cellRange, int cell, int bamFileIndex, char flag, String uniqueId){
			
			for(int i = 0; i < junctions.length; i++){
				if(junctions[i].getStart() == start && junctions[i].getEnd() + 1 == end && (uniqueId.equals("n") || (uniqueId.equals("u") && !junctions[i].getCallsID().contains(ID)))){
					if(flag == 's' || flag == 'l')
						junctions[i].getCells().get(cellRange)[cell]++;
					else
						junctions[i].getBamFiles()[bamFileIndex]++;
					junctions[i].getCallsID().put(ID, 0);
					return true;
				}
			}
			return false;
		
		}
		
		//add 1 exonCall to specific cell
		//cellRange and cell - for singleCell, bamFieIndex - for population
		public void addToExon(int start, int end, String ID, int cellRange, int cell, int index, int bamFileIndex, char flag, String uniqueId){
			if(uniqueId.equals("n") || (uniqueId.equals("u") && !exons[index].getCallsID().contains(ID))){
				if(flag =='s')
					exons[index].getCells().get(cellRange)[cell]++;
				else
					exons[index].getBamFiles()[bamFileIndex]++;
				exons[index].getCallsID().put(ID, 0);

			
			}	
		}
		
		//add 1 intronCall to specific cell
		//cellRange and cell - for singleCell, bamFieIndex - for population
		public void addToIntron(int start, int end, String ID, int cellRange, int cell, int index, int bamFileIndex, char flag, String uniqueId){
			if(uniqueId.equals("n") || (uniqueId.equals("u") && !introns[index].getCallsID().contains(ID))){
				if(flag == 's' || flag == 'l')
					introns[index].getCells().get(cellRange)[cell]++;
				else
					introns[index].getBamFiles()[bamFileIndex]++;
				introns[index].getCallsID().put(ID, 0);
			
			}	
		}
		
		
		//create array of exons
		public void setExonsArray() {
			this.exonsArray = new int[exons.length][2];
			for(int i = 0; i < exons.length; i++){
				exonsArray[i][0] = exons[i].getStart();
				exonsArray[i][1] = exons[i].getEnd();
			}
		}

		
		//create array of introns
		public void setIntronsArray() {
			this.intronsArray = new int[introns.length][2];
			for(int i = 0; i < introns.length; i++){
				intronsArray[i][0] = introns[i].getStart();
				intronsArray[i][1] = introns[i].getEnd();
			}
		}
	//TODO: change the function for hashtable
	public void sortJunctionsByCounter(){
	/*int n=junctions.size();
	LinkedList<Junction> sortedList=new LinkedList<Junction>();
	int[] counterArray=new int[n];
	Junction[] temp=new Junction[n];
	for (int i=0; i<n; i++){
		if(junctions.get(i).getStartCounter()!=0)
			counterArray[junctions.get(i).getStartCounter()-1]++;
		else
			counterArray[junctions.get(i).getEndCounter()-1]++;
	}
	for(int i=1; i<n; i++){
		counterArray[i]=counterArray[i]+counterArray[i-1];
		
	}
	for(int i=n-1; i>=0; i--){
		if(junctions.get(i).getStartCounter()!=0){
			temp[counterArray[junctions.get(i).getStartCounter()-1]-1]=new Junction(junctions.get(i));
			counterArray[junctions.get(i).getStartCounter()-1]=counterArray[junctions.get(i).getStartCounter()-1]-1;
		}
		else{
			temp[counterArray[junctions.get(i).getEndCounter()-1]-1]=new Junction(junctions.get(i));
			counterArray[junctions.get(i).getEndCounter()-1]=counterArray[junctions.get(i).getEndCounter()-1]-1;
		}
	}
	
	for(int i=0; i<temp.length; i++){
		if(temp[i]!=null)
			sortedList.add(temp[i]);
	}
	junctions=sortedList;
	*/
		
	}

	public int sumUnknownJunctionCalls(){
		int sum=0;
		for(int i = 0; i < unknowncells.size(); i++)
			for(int j = 0; j < unknowncells.get(i).length; j++)
				sum += unknowncells.get(i)[j];
		return sum;
	}
	
	public int sumDifCigarCalls(){
		int sum = 0;
		for(int i = 0; i < difCigarCells.size(); i++)
			for(int j = 0; j < difCigarCells.get(i).length; j++)
				sum += difCigarCells.get(i)[j];
		return sum;
	}

	//update the cells counter according to one Call
	public void updateCallSingleCell(String[] line, List<cellsGroup> cells, int startingPosition, int endingPosition, String uniqueId){
		
		Call call = new Call(line);
		String cigar = call.getCigar();

		int cell = call.getCell();
		
		int start = call.getStart();
		String ID = call.getID();

		int cellRange = cellsGroup.findCellsRange(cells, call.getCell());
		if(cellRange == -1){
			System.out.println("Call " + ID + " falling on invalid cell");
			return;
		}
		cell = cell - cells.get(cellRange).getFirst();

		this.cells.get(cellRange)[cell]++;	
		if(!cigar.equals("*")){
			int indexN = findIndex(cigar); //index of first junction
			int indexM = cigar.indexOf((int)'M'); //M index in the cigar
			if((indexM < indexN && indexN != -1) || indexN == -1)
				startWithM(cigar, cell, start, ID, startingPosition, endingPosition, cellRange, -1, 's', uniqueId);
			else{
				int length = Util.convertToInteger(cigar.substring(0, indexN));
				callOnJunction(start, length, ID, cell, cellRange, -1, 's', uniqueId);		
				cigar = cigar.substring(indexN+1);
				startWithM(cigar, cell, start, ID, startingPosition, endingPosition, cellRange, -1, 's', uniqueId);	
			}
							
		}
		else
			this.difCigarCells.get(cellRange)[cell]++;//count the calls with cigar "*"
	}
	
	
	private void startWithM(String cigar, int cell, int start, String ID,int startingPosition, int endingPosition, int cellRange, int bamFileIndex, char flag, String uniqueId){
		
		try {
			int indexM = cigar.indexOf((int) 'M'); //M index in the cigar
			int length = Util.convertToInteger(cigar.substring(0, indexM)); //length on the call until the M


			//if the call starts before the beginning of the gene
			String[] temp = startBeforeGene(start, startingPosition, length, ID, cell, cellRange, indexM, cigar, bamFileIndex, flag, uniqueId);
			cigar = temp[0];
			start = Util.convertToInteger(temp[1]);
			indexM = Util.convertToInteger(temp[2]);

			//create introns and exons array
			int[][] intronsArray = getIntronsArray();
			int[][] exonsArray = getExonsArray();

			//initial first indexes
			int exonIndex = 0, intronIndex = 0, firstExonAdded, firstIntronAdded, exonEnd = 0, intronEnd = 0, intronStart = 0;
			int exonStart = exonsArray[exonIndex][0];
			if (intronsArray.length > 0)
				intronStart = intronsArray[intronIndex][0];

			while (!cigar.equals("")) {//while the call did not finished
				firstExonAdded = -1;
				firstIntronAdded = -1;

				//check on which exons the call falls
				while (exonIndex < exonsArray.length && (exonStart <= start || exonStart == exonsArray[exonIndex][0] || (startingPosition > start))) {
					exonStart = exonsArray[exonIndex][0];
					exonEnd = exonsArray[exonIndex][1];

					if ((start <= exonEnd && start >= exonStart) || (start + length <= exonEnd && start + length > exonStart)) {
						addToExon(exonStart, exonEnd, ID, cellRange, cell, exonIndex, bamFileIndex, flag, uniqueId);

						if (firstExonAdded == -1)
							firstExonAdded = exonIndex;
					}
					exonIndex++;
				}
				if (exonIndex > 0)
					exonIndex--;

				//check on which introns the call falls
				while (intronIndex < intronsArray.length && (intronStart <= start || intronStart == intronsArray[intronIndex][0])) {
					intronStart = intronsArray[intronIndex][0];
					intronEnd = intronsArray[intronIndex][1];
					if ((start <= intronEnd && start >= intronStart) || (start + length <= intronEnd && start + length > intronStart)) {
						addToIntron(intronStart, intronEnd, ID, cellRange, cell, intronIndex, bamFileIndex, flag, uniqueId);
						if (firstIntronAdded == -1)
							firstIntronAdded = intronIndex;
					}
					intronIndex++;
				}
				if (intronIndex > 0)
					intronIndex--;


				if (firstExonAdded != -1)
					exonIndex = firstExonAdded;
				if (firstIntronAdded != -1)
					intronIndex = firstIntronAdded;


				//check on which junction the call falls
				if (exonIndex == exonsArray.length)
					break;
				//count the junction
				start = start + length;
				if (cigar.length() > indexM) {
					cigar = cigar.substring(indexM + 1);
					int indexN = findIndex(cigar); //N index in the cigar
					if (indexN == -1)
						break;
					length = Util.convertToInteger(cigar.substring(0, indexN));
					if (length > 3) {
						callOnJunction(exonsArray[exonIndex][1] + 1, length, ID, cell, cellRange, bamFileIndex, flag, uniqueId);
						start = start + length;
						if (cigar.length() == indexN)
							break;
						cigar = cigar.substring(indexN + 1);
						indexM = cigar.indexOf((int) 'M');
						if (indexM == -1)
							break;
						length = Util.convertToInteger(cigar.substring(0, indexM));
					} else {
						indexM = cigar.indexOf((int) 'M');
						if (indexM == -1)
							break;
						length = Util.convertToInteger(cigar.substring(0, indexN)) + Util.convertToInteger(cigar.substring(indexN + 1, indexM));

					}


				}
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}

		
	}
	
	private void callOnJunction(int start, int length,String ID, int cell, int cellRange, int bamFileIndex, char flag, String uniqueId){		
		boolean found = addToJunction(start, start+length, ID, cellRange, cell, bamFileIndex, flag, uniqueId);
		if(!found && flag == 's')
			this.getUnknowncells().get(cellRange)[cell]++;
		if(!found && flag == 'p')
			this.unknownBF[bamFileIndex]++;

	}
	
	private String[] startBeforeGene(int start, int startingPosition, int length, String ID, int cell, int cellRange, int indexM, String cigar, int bamFileIndex, char flag, String uniqueId){
		int lengthBeforStart = startingPosition-start;
		while(start < startingPosition){
			if(length > lengthBeforStart)
				addToExon(getExonsArray()[0][0], getExonsArray()[0][1], ID, cellRange, cell, 0, bamFileIndex, flag, uniqueId);
			start = start + length;
			lengthBeforStart = lengthBeforStart-length;
			if(cigar.length() == indexM)
				break;
			else
				cigar = cigar.substring(indexM + 1);
			int indexN = findIndex(cigar);
			if(indexN == -1)
				break;
			length = Util.convertToInteger(cigar.substring(0, indexN));
			if(lengthBeforStart > length){
				start = start + length;
				lengthBeforStart = lengthBeforStart - length;
			}			
			else
				start = start + length;
				
			
			if(cigar.length() == indexN)
				break;
			else
				cigar = cigar.substring(indexN + 1);
			indexM = cigar.indexOf((int)'M');
			if(indexM == -1)
				break;			
			}
		String[] returnArr = {cigar, Integer.toString(start), Integer.toString(indexM)};
		return returnArr;
	}
	
	//update the cells counter according to one Call
	public void updateCallPopulation(String[] line, int startingPosition, int endingPosition, int BamIndex, String uniqueId){
		Call call = new Call(line);
		String cigar = call.getCigar();
		int start = call.getStart();
		String ID = call.getID();
		totalCallsBF[BamIndex]++;
		if(!cigar.equals("*")){
			int indexN = findIndex(cigar); //index of first junction
			int indexM = cigar.indexOf((int)'M'); //M index in the cigar
			if((indexM < indexN && indexN != -1) || indexN == -1)
				startWithM(cigar,-1, start, ID, startingPosition, endingPosition,-1, BamIndex,'p', uniqueId);//-1 is for singleCell arguments
			else{
				int length = Util.convertToInteger(cigar.substring(0, indexN));
				callOnJunction(start, length, ID,-1, -1, BamIndex, 'p', uniqueId);//-11 i for singleCell arguments		
				cigar = cigar.substring(indexN+1);
				startWithM(cigar,-1, start, ID, startingPosition, endingPosition,-1, BamIndex, 'p', uniqueId);//-1 is for singleCell arguments
			}
							
		}
		else 
			this.difCigarBF[BamIndex]++;
		
	}

	public void writeGeneToFile(BufferedWriter file){
		try {
			// Write gene name and chrom
			file.write(name + '\t' + chrom + '\t');

			// Write cells
			if(cells != null) {
				for (int i = 0; i < cells.size(); i++) {
					file.write("1-" + cells.get(i).length);
					if (i < cells.size() - 1)
						file.write(',');
				}
				file.write('\t');
			}
			else
				file.write("-" + '\t');

			// Write exon starts
			for(int i = 0; i < exonStarts.size(); i++){
				file.write(exonStarts.get(i));
				if(i < exonStarts.size() - 1)
					file.write('/');
			}
			file.write('\t');

			// Write exon ends
			for(int i = 0; i < exonEnds.size(); i++){
				file.write(exonEnds.get(i));
				if(i < exonEnds.size() - 1)
					file.write('/');
			}
			file.write('\t');
			file.write(String.valueOf(totalCallsBF.length));
		}
		catch (IOException e){
			e.printStackTrace();
		}
	}

	public static Gene createGeneFromFile(String geneLineInFile){
		String[] splitLine = geneLineInFile.split("\t");
		String name = splitLine[0];
		String chrom = splitLine[1];

		//Create cells list
		List<cellsGroup> cells = new ArrayList<cellsGroup>();
		if(!splitLine[2].equals("-")) {
			String[] cellsStrings = splitLine[2].split(",");
			for (int i = 0; i < cellsStrings.length; i++) {
				String[] cellStartAndEnd = cellsStrings[i].split("-");
				cellsGroup cell = new cellsGroup(Integer.parseInt(cellStartAndEnd[0]), Integer.parseInt(cellStartAndEnd[1]));
				cells.add(cell);
			}
		}

		//Create exons starts list
		LinkedList<String> starts = new LinkedList<>();
		String[] exonStarts = splitLine[3].split("/");
		for(String start: exonStarts){
			starts.add(start);
		}
		//Create exons ends list
		LinkedList<String> ends = new LinkedList<>();
		String[] exonsEnds = splitLine[4].split("/");
		for(String end: exonsEnds){
			ends.add(end);
		}

		int bamFilesNumber = Integer.parseInt(splitLine[5]);

		Gene gene = new Gene(name, chrom, cells, starts, ends, bamFilesNumber);
		return gene;
	}
}
