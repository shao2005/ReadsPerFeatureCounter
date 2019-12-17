import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.*;
import java.io.FileWriter;
import java.util.concurrent.TimeUnit;

import org.omg.CORBA.portable.InputStream;
import org.omg.CORBA.portable.OutputStream;

public class Processing {

	/**
	 * main function
	 * 
	 * @param cells
	 *            List of cells ranges (SingleCell)
	 * @param BAMFile
	 *            bamFile name (singleCell)
	 * @param BamFolder
	 *            folder name of the bam file folder (singleCell)
	 * @param genesFile
	 *            transcripts file
	 * @param outputFolder
	 *            folder name for all the output files
	 * @param flag
	 *           'p' if population, 's' if singleCell
	 * @param bamFiles
	 *            array of bamFiles names (population)
	 * @param pathForPopulationFolder
	 *            path in order to get to the bam files int the array
	 *            (population)
	 */
	public void process(List<cellsGroup> cells, String BAMFile, String BamFolder, String genesFile,
			String outputFolder, char flag, String[] bamFiles, String pathForPopulationFolder, String geneToStart,
			String UniqueId, String picardPath) {
		IO io = new IO();

		try {
			//create output folder
			if((flag == 's' || flag == 'p') && geneToStart.equals("-")){
				File theDir = new File(outputFolder);
				File LG = new File(outputFolder + "/LG");
				theDir.mkdir();
				LG.mkdir();
			}
			
			FileWriter runLater = io.createRunLaterFile(outputFolder, "/runLater.txt");
			boolean append = false;
			if(!geneToStart.equals("-"))
				append = true;
			
			// files creation
			FileWriter trenscriptsFile = new FileWriter(outputFolder + "/junctions.csv", append);
			FileWriter exonsFile = new FileWriter(outputFolder + "/exons.csv", append);
			FileWriter intronsFile = new FileWriter(outputFolder + "/introns.csv", append);
			FileWriter unknownFile = new FileWriter(outputFolder + "/unknowenJunctions.csv", append);
			FileWriter geneCallsFile = new FileWriter(outputFolder + "/geneCalls.csv", append);
			FileWriter difCigarFile = new FileWriter(outputFolder + "/difCigar.csv", append);
			FileWriter times = new FileWriter(outputFolder + "/timestamps.txt", append);
			FileWriter[] files = { trenscriptsFile, exonsFile, intronsFile, unknownFile, geneCallsFile, difCigarFile, times};

			io.createJunctionsFile(trenscriptsFile, outputFolder + "/junctions.csv", cells, bamFiles, flag);
			io.createExonsFile(exonsFile, outputFolder + "/exons.csv", cells, bamFiles, flag);
			io.createIntronsFile(intronsFile, outputFolder + "/introns.csv", cells, bamFiles, flag);
			io.createUnknownJunctionFile(unknownFile, outputFolder + "/unknowenJunctions.csv", cells, bamFiles, flag);
			io.createGeneCallsFile(geneCallsFile, outputFolder + "/geneCalls.csv", cells, bamFiles, flag);
			io.createUnknownJunctionFile(difCigarFile, outputFolder + "/difCigar.csv", cells, bamFiles, flag);
			io.closeFiles(files);
			BufferedReader br = new BufferedReader(new FileReader(genesFile));
			// end files creation
			
			HashMap<String, Gene> genes = io.readTranscriptsFileIntoHashMap(br, cells, bamFiles.length);
			io.writeGenesHashSetIntoFile(outputFolder + "/genesObjects.txt", genes);
			int index = 0;
			
			List<String> needToRun = null;
			if(flag == 'l' || flag == 'L')
				needToRun = fileToList(outputFolder.substring(0, outputFolder.length()-3) + "/runLater.txt");


			BufferedReader genesObjectsFile = new BufferedReader(new FileReader(outputFolder + "/genesObjects.txt"));
			boolean foundFirstGene = false;
			String geneLine = "";
			while((geneLine = genesObjectsFile.readLine()) != null){
				Gene gene = Gene.createGeneFromFile(geneLine);
				String geneName = gene.getName();
				index++;
				//if program need to start from specific gene, skip genes until get to the gene
				if(!foundFirstGene && !geneToStart.equals("-")){
					if(geneName.equals(geneToStart))
						foundFirstGene = true;
					else
						continue;
				}
				
				if (flag == 's' || flag == 'p' || ((flag == 'l' || flag == 'L') && needToRun.contains(geneName))) {
					// open files to write
					trenscriptsFile = new FileWriter(outputFolder + "/junctions.csv", true);
					exonsFile = new FileWriter(outputFolder + "/exons.csv", true);
					intronsFile = new FileWriter(outputFolder + "/introns.csv", true);
					unknownFile = new FileWriter(outputFolder + "/unknowenJunctions.csv", true);
					geneCallsFile = new FileWriter(outputFolder + "/geneCalls.csv", true);
					difCigarFile = new FileWriter(outputFolder + "/difCigar.csv", true);
					times = new FileWriter(outputFolder + "/timestamps.txt", true);
					runLater = new FileWriter(outputFolder + "/runLater.txt", true);
	
					FileWriter[] files2 = { trenscriptsFile, exonsFile, intronsFile, unknownFile, geneCallsFile,
							difCigarFile, times, runLater };
	
					System.out.println(index + "." + geneName + " ");
					// write the gene to timestamps file
					times.write(index + ". " + gene.getName() + "	");
					times.write(new Timestamp(new java.util.Date().getTime()).toString());
					times.append(System.lineSeparator());
					List<List<String>> starts = createStartsList(gene);
					List<List<String>> ends = createEndsList(gene);
					if (starts != null) {
	
						gene.setExons(starts, ends, cells, bamFiles.length, flag);
						gene.setIntrons(starts, ends, cells,bamFiles.length, flag);
						gene.setJunctions(cells, bamFiles.length, flag);
						gene.setExonsArray();
						gene.setIntronsArray();
						int startingPosition = Util.findMin(starts);
						int endingPosition = Util.findMax(ends);
						//changed 5th parameter from BAMFile to bamFiles[0]
						io.createIntervalList(gene, startingPosition, endingPosition, pathForPopulationFolder, BAMFile, outputFolder);
						int callsNumber = 0;
						
						if (flag == 's') {
							//changed 1st parameter from BAMFile to bamFiles[0]
							readInterval(BAMFile, outputFolder, pathForPopulationFolder, picardPath, flag);
							callsNumber = updateCells(BAMFile, gene, cells, startingPosition, endingPosition,runLater, flag, -1, outputFolder, UniqueId); // the one before last argument is connected to population code
						}
						else if(flag == 'l'){
							//changed 1st parameter from BAMFile to bamFiles[0]
							readInterval(BAMFile, outputFolder, pathForPopulationFolder, picardPath, flag);
							callsNumber = updateCellsForLongGenes(BAMFile, outputFolder, gene, cells, startingPosition,endingPosition, 's', 0, UniqueId); // the last argument is connected to population code
						}
						else {
							for (int j = 0; j < bamFiles.length; j++) {
								readInterval(bamFiles[j], outputFolder, pathForPopulationFolder, picardPath, flag);
								if(flag == 'p')
									callsNumber = updateCells(bamFiles[j], gene, cells, startingPosition, endingPosition, runLater, flag, j, outputFolder, UniqueId);
								else
									callsNumber = updateCellsForLongGenes(bamFiles[j], outputFolder, gene, cells, startingPosition,endingPosition, 'p', j, UniqueId);
								if(callsNumber == -1)
									break;										
							}
						}

						if(callsNumber != -1){//write to file only if gene finished
							io.writeJunction(trenscriptsFile, gene, cells, callsNumber, flag);
							io.writeExon(exonsFile, gene, cells, callsNumber, flag);
							io.writeIntron(intronsFile, gene, cells, callsNumber, flag);
							io.writeUnknownJunction(unknownFile, gene, cells, startingPosition, endingPosition, callsNumber, flag);						
							io.writeGeneCell(geneCallsFile, gene, cells, startingPosition, endingPosition, callsNumber, callsNumber, flag);										
							io.writedifCigar(difCigarFile, gene, cells, startingPosition, endingPosition, callsNumber, flag);
						}
					}
					io.closeFiles(files2);
				}
			}
			trenscriptsFile.close();
			exonsFile.close();
			intronsFile.close();
			unknownFile.close();
			difCigarFile.close();
			geneCallsFile.close();
			times.close();
			br.close();

		} catch (IOException e) {
			e.printStackTrace();
			System.out.println('e');
		}
	}
	
	
	// reads the BAM file and copy to 'one_interval' file the calls that fall on
	// the gene interval

	public boolean readInterval(String file, String outputFolder, String pathToBamFolder, String picardPath, char flag) {
		try {// running on windows
			if(!new File(pathToBamFolder, file).exists())
				System.out.println("The bam file " + file + " does not exists inside " + pathToBamFolder);
			if(!new File(picardPath, "picard.jar").exists())
				System.out.println("The file picard.jar does no exists inside " + picardPath);
			if (System.getProperty("os.name").compareTo("Windows 7") == 0 || System.getProperty("os.name").compareTo("Windows 10") == 0) {
				String command = "cmd /c start /wait cmd.exe /K \"java -jar " + picardPath + "\\picard.jar ViewSam I="+ pathToBamFolder +"\\"+ file + " INTERVAL_LIST=" + outputFolder + "\\one_gene_";
				if(flag == 'l' || flag == 's')		
					command += file + ".interval_list > " + outputFolder + "\\one_gene_" + file + ".txt";
				else
					command += ".interval_list > " + outputFolder + "\\one_gene_" + file + ".txt";

				Process p = Runtime.getRuntime().exec(command);
				p.waitFor();
				Runtime.getRuntime().exec("taskkill /F /IM cmd.exe");
				//Thread.sleep(1000);
			} else {// running on cluster
				String command = "java -jar " + picardPath + "/picard.jar ViewSam I="+ pathToBamFolder +"/"+ file + " INTERVAL_LIST=" + outputFolder + "/one_gene_";
				if(flag == 'l' || flag == 's' || flag =='p')
					command += file + ".interval_list > " + outputFolder + "/one_gene_" + file + ".txt";
				else
					command += ".interval_list > " + outputFolder + "/one_gene_" + file + ".txt";
				String[] cmd = { "/bin/sh", "-c", command };
				Process p = Runtime.getRuntime().exec(cmd);
				p.waitFor();
                Runtime.getRuntime().exec("kill $p");
			}
		} catch (Exception e) {
			System.out.println(e);
		}
		return true;
	}

	// update the cells counter according to the calls
	public int updateCells(String fileName, Gene gene, List<cellsGroup> cells, int startingPosition, int endingPosition,FileWriter checkLater, char flag, int bamIndex, String outputFolder, String uniqueId) {
		String line = null;
		int countLines = 0;
		int checkIfFileIsEmpty=0;
		try {

			BufferedReader br = new BufferedReader(new FileReader( outputFolder + "/one_gene_" + fileName + ".txt"));
			long startTime = System.currentTimeMillis();
			while ((line = br.readLine()) != null) {
				checkIfFileIsEmpty++;
				String[] call = line.split("\t");
				if (call[0].substring(0, 1).equals("@") || call.length == 1) // check if we need to ignore from this line
					continue;
				else {
                    countLines++;
					if (flag == 's') 
						gene.updateCallSingleCell(call, cells, startingPosition, endingPosition, uniqueId);
					
					else
						gene.updateCallPopulation(call, startingPosition, endingPosition, bamIndex, uniqueId);
					
					long currTime = System.currentTimeMillis();
					if (currTime - startTime > 120000) { // check if the gene did not finish in time

						checkLater.write(gene.getName());
						checkLater.append(System.lineSeparator());
						return -1;
					}
				}
			}
			if(checkIfFileIsEmpty == 0){
				System.out.println("The calls file is empty. Please check your parameters.");
				return -2;
			}
			br.close();
		}

		catch (IOException e) {
			e.printStackTrace();
		}
		return countLines;
	}

	// update the cells counter according to the calls. no need for flag because
	public int updateCellsForLongGenes(String fileName, String outputFolder, Gene gene, List<cellsGroup> cells, int startingPosition,int endingPosition, char flag, int bamIndex, String uniqueId) {
		String line = null;
		int countLines = 0;
		int checkIfFileIsEmpty=0;

		try {
			BufferedReader br = new BufferedReader(new FileReader(outputFolder + "/one_gene_" + fileName + ".txt"));
			while ((line = br.readLine()) != null) {
				checkIfFileIsEmpty++;
				String[] call = line.split("\t");
				if (call[0].substring(0, 1).equals("@")) // check if we need to ignore from this line
					continue;
				else {
					countLines++;
					if (flag == 's') 
						gene.updateCallSingleCell(call, cells, startingPosition, endingPosition, uniqueId);
					
					else
						gene.updateCallPopulation(call, startingPosition, endingPosition, bamIndex, uniqueId);
				}
			}
			if(checkIfFileIsEmpty == 0){
				System.out.println("The calls file is empty. Please check your parameters.");
				return -2;
			}
			br.close();
		}

		catch (IOException e) {
			e.printStackTrace();
		}
		return countLines;
	}

	// insert genes from file to list
	public List<String> fileToList(String fileName) {
		List<String> genesList = new ArrayList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line = "";
			while ((line = br.readLine()) != null)
				genesList.add(line);

		} catch (Exception e) {
			e.printStackTrace();
		}
		return genesList;

	}

	// create a list of the exons starts
	public List<List<String>> createStartsList(Gene gene) {
		List<List<String>> starts = new ArrayList<List<String>>();
		for (int j = 0; j < gene.getExonStarts().size(); j++) {
			String[] start = gene.getExonStarts().get(j).split(",");
			starts.add(Arrays.asList(start));
		}
		if (gene.getName().equals("NCR2")) {
			String[] start = new String[] { "41303357", "41303824", "41309531", "41309773", "41310655" };
			starts.add(Arrays.asList(start));
		}

		return starts;
	}

	// create a list of the exons ends
	public List<List<String>> createEndsList(Gene gene) {
		List<List<String>> ends = new ArrayList<List<String>>();
		for (int j = 0; j < gene.getExonEnds().size(); j++) {
			String[] end = gene.getExonEnds().get(j).split(",");
			ends.add(Arrays.asList(end));
		}
		if (gene.getName().equals("NCR2")) {
			String[] end = new String[] { "41303666", "41304166", "41309667", "41309887", "41311139" };
			ends.add(Arrays.asList(end));
		}

		return ends;
	}

	// combine between the 2 files and find the gene name according to the
	// symbol. if cound'nt find, return null.
	public String findName(LinkedList<Pair> names, String symbol) {
		for (int i = 0; i < names.size(); i++)
			if (names.get(i).getSign().equals(symbol))
				return names.get(i).getName();
		return null;

	}

	// count lines in file
	public int countLinesInFile(String filename) throws IOException {
		BufferedInputStream is = new BufferedInputStream(new FileInputStream(filename));
		try {
			byte[] c = new byte[1024];
			int count = 0;
			int readChars = 0;
			boolean endsWithoutNewLine = false;
			while ((readChars = is.read(c)) != -1) {
				for (int i = 0; i < readChars; ++i) {
					if (c[i] == '\n')
						++count;
				}
				endsWithoutNewLine = (c[readChars - 1] != '\n');
			}
			if (endsWithoutNewLine) {
				++count;
			}
			return count;
		} finally {
			is.close();
		}
	}

	// create a list of pairs of the gene name and sign
	public LinkedList<Pair> createPairsList1(String file, int nameIndex) {
		LinkedList<Pair> tempNames = new LinkedList<>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] fullLine = line.split("\t");
				Pair pair = new Pair(fullLine[0], fullLine[4]);
				tempNames.add(pair);
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return tempNames;
	}

	public Map<String, String> createPairsMap(String file, int nameIndex) {
		Map<String, String> m = new HashMap<String, String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] fullLine = line.split("\t");
				m.put(fullLine[0], fullLine[nameIndex]);
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return m;
	}
	
	public static void readGTF(String path, String gtfFile){
		try {
			BufferedReader br = new BufferedReader(new FileReader(gtfFile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(path + "gtfFile.txt"));
			String line;

			// ignore first 4 lines
			for (int i = 0; i < 5; i++) 
				line = br.readLine();
			int index=0;
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				String[] transcript = row[8].split(";");
				if(transcript.length > 4){
					String name = transcript[4].substring(12, transcript[4].length() - 1);
					if(name.equals("NIPAL3"))
						bw.write(index + ". " + line + "\n");
				}
				index++;

			}
			br.close();
			bw.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}
			
	}

	// create transcripts list according to "Homo_sapiens.Ensembl.GRCh38.82.gtf"
		// transcripts file
		public Hashtable<String, Gene> createGenesListHuman38(String gtfFile) {
			Hashtable<String, Gene> genes = new Hashtable<String, Gene>();
			Hashtable<String, String> idToName = new Hashtable<String, String>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(gtfFile));
				String line;

				// ignore first 4 lines
				for (int i = 0; i < 5; i++) {
					line = br.readLine();
				}

				String[] row;
				int index=0;
				while ((line = br.readLine()) != null) {
					row = line.split("\t");
					String name = null, chrom, id = null;
					chrom = row[0];
					
					if (!row[2].equals("gene")) {
						row = line.split("\t");
						while (row[2].equals("transcript") && line != null) {
							String exonStart = "", exonEnd = "";
							String[] transcript = row[8].split(";");
							//look for gene name
							for(int i = 1; i < transcript.length && name == null; i++){
								if(transcript[i].contains("gene_name"))
									name = transcript[i].substring(12, transcript[4].length() - 1);
							}

							
							id = transcript[0].substring(9, transcript[0].length() - 1);
							while ((line = br.readLine()) != null) {
								String[] exon = line.split("\t");
								if (exon[2].equals("exon")) {
									exonStart = exonStart + exon[3] + ",";
									exonEnd = exonEnd + exon[4] + ",";
								} else {
									row = exon;
									break;
								}

							}
							if (!idToName.containsKey(id)) {
								System.out.println("{" +index+"."+ name + "}");
								idToName.put(id, name);
								index++;
								LinkedList<String> exonStarts = new LinkedList<String>();
								LinkedList<String> exonEnds = new LinkedList<String>();
								//exonStarts.add(exonStart);
								//exonEnds.add(exonEnd);
								genes.put(name, new Gene(name, chrom, null, exonStarts, exonEnds, 0));

							}
							name = idToName.get(id);
							genes.get(name).getExonStarts().add(exonStart);
							genes.get(name).getExonEnds().add(exonEnd);

						}

					}

				}
				br.close();

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}

			return genes;

		}
		
		public Hashtable<String, Gene> createGenesListHuman(String path, String gtfFile){
			Hashtable<String, Gene> genes = new Hashtable<>();
			Hashtable<String, String> idToName = new Hashtable<>();
			Hashtable<String, String> geneChrom = new Hashtable<>();
			Hashtable<String, Set<String>> geneTrancripts = new Hashtable<>();
			Hashtable<String, List<String>> transcriptStarts = new Hashtable<>();
			Hashtable<String, List<String>> transcriptEnds = new Hashtable<>();
			HashSet<String> genesToDelete = new HashSet<>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(gtfFile));
				BufferedWriter deletedGenesFile = new BufferedWriter(new FileWriter(path + "doubleGenes.txt"));
				String line;
				
				// ignore first 4 lines
				for (int i = 0; i < 5; i++) 
					line = br.readLine();
				
				
				line = br.readLine();
				while(line != null){
					String[] row = line.split("\t");
					String chrom = row[0];
					String geneSign = row[6];
					
					//search for genes lines 
					if(row[2].equals("gene")){
						String[] rowDetails = row[8].split(";");
						
						//initial gene						
						String geneName = rowDetails[2].substring(12, rowDetails[2].length() - 1);
						String geneId = rowDetails[0].substring(9, rowDetails[0].length() - 1);
						
						System.out.println(geneName + " ");
						if(!geneChrom.containsKey(geneName) && !genesToDelete.contains(geneName)){		
							geneChrom.put(geneName, chrom);
							geneTrancripts.put(geneName, new HashSet<String>());
							idToName.put(geneId, geneName);
						}//check if the gene has another chromosome
						else if(geneChrom.containsKey(geneName) && !geneChrom.get(geneName).equals(chrom)){
							genesToDelete.add(geneName);
							geneChrom.remove(geneName);
						}
						
						line = br.readLine();
												
						while(line != null){
							row = line.split("\t");
							//check only exon lines
							if(row[2].equals("exon")){
								String transcriptSign = row[6];
								rowDetails = row[8].split(";");
								String geneIdForExon = rowDetails[0].substring(9, rowDetails[0].length() - 1);
								String transcriptID = null;
								for(int i = 0; i < rowDetails.length && transcriptID == null; i++)
									if(rowDetails[i].contains("transcript_id"))
										transcriptID = rowDetails[i].substring(16, rowDetails[i].length() - 1);
								
							
								if(!transcriptStarts.containsKey(transcriptID)){
									transcriptStarts.put(transcriptID, new ArrayList<>());
									transcriptEnds.put(transcriptID, new ArrayList<>());
								}
								
								//add exon to transcript
								if(transcriptSign.equals(geneSign)){
									geneTrancripts.get(geneName).add(transcriptID);
									if(!idToName.containsKey(geneIdForExon))
										idToName.put(geneIdForExon, geneName);
								}
								else if(idToName.containsKey(geneIdForExon))
									geneTrancripts.get(idToName.get(geneIdForExon)).add(transcriptID);
								
								transcriptStarts.get(transcriptID).add(row[3]);
								transcriptEnds.get(transcriptID).add(row[4]);

															
								line = br.readLine();
																			
							}
							else if(row[2].equals("transcript") || row[2].equals("CDS") || row[2].equals("start_codon") || row[2].equals("stop_codon") || row[2].equals("five_prime_utr") || row[2].equals("three_prime_utr") || row[2].equals("Selenocysteine"))
								line = br.readLine();
							else if(row[2].equals("gene"))
								break;
							else
								System.out.println(line);

						}
					}
				}
				List<String> transToRemove = new ArrayList<>();

				//create genes according to transcripts
				for(String name : geneChrom.keySet()){
					String chrom = geneChrom.get(name);
					LinkedList<String> exonStarts = new LinkedList<String>();
					LinkedList<String> exonEnds = new LinkedList<String>();	
					for(String transcript : geneTrancripts.get(name)){
						exonStarts.add(ListToString(transcriptStarts.get(transcript)));
						exonEnds.add(ListToString(transcriptEnds.get(transcript)));	
					}
		
					genes.put(name, new Gene(name, chrom, null, exonStarts, exonEnds, 0));
						
				}
				for(int i = 0; i < transToRemove.size(); i++)
					transcriptStarts.remove(transToRemove.get(i));
				
				//write all deleted genes to file
				for(String gene: genesToDelete)
					deletedGenesFile.write(gene + "\n");

				
				deletedGenesFile.close();
				br.close();
					
			}
			catch(IOException e){
				e.printStackTrace();
			}
				
			return genes;
		}
		
		private String ListToString(List<String> list){
			String ans = "";
			for(int i = 0; i < list.size(); i++)
				ans += list.get(i) + ",";
			return ans;
		}
		
		// create transcripts list according to
		// "mm10_combined_annotated_oneName.gtf" transcripts file
 		public Hashtable<String, Gene> createGenesListMouse(String path, String gtfFile) {
			Hashtable<String, Gene> genes = new Hashtable<String, Gene>();
			Hashtable<String, String> geneChrom = new Hashtable<>();
			Hashtable<String, Set<String>> geneTrancripts = new Hashtable<>();
			Hashtable<String, List<String>> transcriptStarts = new Hashtable<>();
			Hashtable<String, List<String>> transcriptEnds = new Hashtable<>();
			HashSet<String> genesToDelete = new HashSet<>();

			try {
				BufferedReader br = new BufferedReader(new FileReader(gtfFile));
				BufferedWriter deletedGenesFile = new BufferedWriter(new FileWriter(path + "doubleGenes.txt"));
				String line;
				
				while((line = br.readLine()) != null){
					String[] row = line.split("\t");
					
					//choose only exon lines
					if (row[2].equals("exon")) {
						String chrom = row[0];
						String[] transcript = row[8].split(";");
						String name = transcript[0].substring(9, transcript[0].length() - 1);
						String transcriptID = transcript[1].substring(16, transcript[1].length() - 1);
						if(!geneChrom.containsKey(name) && !genesToDelete.contains(name)){
							geneChrom.put(name, chrom);
							geneTrancripts.put(name,  new HashSet<String>());
							System.out.println(name);
						}//check if the gene has another chromosome
						else if(geneChrom.containsKey(name) && !geneChrom.get(name).equals(chrom)){
							genesToDelete.add(name);
							geneChrom.remove(name);
							continue;
						}		
						
						geneTrancripts.get(name).add(transcriptID);
						
						if(!transcriptStarts.containsKey(transcriptID)){
							transcriptStarts.put(transcriptID, new ArrayList<>());
							transcriptEnds.put(transcriptID, new ArrayList<>());
						}
						
						//add exon to transcript
						transcriptStarts.get(transcriptID).add(row[3]);
						transcriptEnds.get(transcriptID).add(row[4]);
					}
				}
				//create genes according to transcripts
				for(String name : geneChrom.keySet()){
					String chrom = geneChrom.get(name);
					LinkedList<String> exonStarts = new LinkedList<String>();
					LinkedList<String> exonEnds = new LinkedList<String>();	
					for(String transcript : geneTrancripts.get(name)){
						exonStarts.add(ListToString(transcriptStarts.get(transcript)));
						exonEnds.add(ListToString(transcriptEnds.get(transcript)));	
					}
		
					genes.put(name, new Gene(name, chrom, null, exonStarts, exonEnds, 0));
						
				}

				
				//write all deleted genes to file
				for(String gene: genesToDelete)
					deletedGenesFile.write(gene + "\n");
				
				deletedGenesFile.close();
				br.close();
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}

			return genes;

		}

		// create genes file
		public void createGenesFile(String destFile, Hashtable<String, Gene> genes, String path) {
			try {
				FileWriter dest = new FileWriter(path + destFile);
				Set<String> keys = genes.keySet();
				for (String key : keys) {
					for (int j = 0; j < genes.get(key).getExonStarts().size(); j++) {
						dest.write(genes.get(key).getName());
						dest.append('\t');
						dest.write(genes.get(key).getChrom());
						dest.append('\t');
						dest.write(genes.get(key).getExonStarts().get(j).substring(0,
								genes.get(key).getExonStarts().get(j).length() - 1));
						dest.append('\t');
						dest.write(genes.get(key).getExonEnds().get(j).substring(0,
								genes.get(key).getExonEnds().get(j).length() - 1));
						dest.append(System.lineSeparator());
					}
				}
				dest.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
}
