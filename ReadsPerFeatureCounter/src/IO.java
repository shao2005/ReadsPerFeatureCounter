import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.io.FileWriter;

public class IO {

	// creates the headers of the junctions file
	public void createJunctionsFile(FileWriter writer, String fileName, List<cellsGroup> cells, String[] BAMFiles,
									char flag) {
		try {
			writer.write("Gene Name");
			writer.append(',');
			writer.write("Chromosome");
			writer.append(',');
			if (flag == 's' || flag == 'l') {
				writer.write("Index");
				writer.append(',');
			}
			writer.write("Starting Point");
			writer.append(',');
			writer.write("Ending Point");
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++)
					for (int j = cells.get(i).getFirst(); j <= cells.get(i).getLast(); j++) {
						writer.append(',');
						writer.write("cell" + j);
					}
			} else {
				for (int i = 0; i < BAMFiles.length; i++) {
					writer.append(',');
					writer.write(BAMFiles[i]);
				}

			}
			if (flag == 's' || flag == 'l') {
				writer.append(',');
				writer.write("Total");
			}
			writer.append('\n');
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// creates the headers of the exons file
	public void createExonsFile(FileWriter writer, String fileName, List<cellsGroup> cells, String[] BAMFiles,
								char flag) {
		try {
			writer.write("Gene Name");
			writer.append(',');
			writer.write("Chromosome");
			writer.append(',');
			writer.write("Exon Start");
			writer.append(',');
			writer.write("Exon End");
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++)
					for (int j = cells.get(i).getFirst(); j <= cells.get(i).getLast(); j++) {
						writer.append(',');
						writer.write("cell" + j);
					}
			} else {
				for (int i = 0; i < BAMFiles.length; i++) {
					writer.append(',');
					writer.write(BAMFiles[i]);
				}

			}
			if (flag == 's' || flag == 'l') {
				writer.append(',');
				writer.write("Total");
			}
			writer.append('\n');
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// creates the headers of the introns file
	public void createIntronsFile(FileWriter writer, String fileName, List<cellsGroup> cells, String[] BAMFiles,
								  char flag) {
		try {
			writer.write("Gene Name");
			writer.append(',');
			writer.write("Chromosome");
			writer.append(',');
			writer.write("Intron Start");
			writer.append(',');
			writer.write("Intron End");
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++)
					for (int j = cells.get(i).getFirst(); j <= cells.get(i).getLast(); j++) {
						writer.append(',');
						writer.write("cell" + j);
					}
			} else {
				for (int i = 0; i < BAMFiles.length; i++) {
					writer.append(',');
					writer.write(BAMFiles[i]);
				}

			}
			if (flag == 's' || flag == 'l') {
				writer.append(',');
				writer.write("Total");
			}
			writer.append('\n');
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// creates the headers of the unknownJunctions file
	public void createUnknownJunctionFile(FileWriter writer, String fileName, List<cellsGroup> cells, String[] BAMFiles,
										  char flag) {
		try {
			writer.write("Gene Name");
			writer.append(',');
			writer.write("Chromosome");
			writer.append(',');
			writer.write("Start");
			writer.append(',');
			writer.write("End");
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++)
					for (int j = cells.get(i).getFirst(); j <= cells.get(i).getLast(); j++) {
						writer.append(',');
						writer.write("cell" + j);
					}
			} else {
				for (int i = 0; i < BAMFiles.length; i++) {
					writer.append(',');
					writer.write(BAMFiles[i]);
				}
			}
			if (flag == 's' || flag == 'l') {
				writer.append(',');
				writer.write("Total");
			}
			writer.append('\n');
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// creates the headers of the geneCalls file
	public void createGeneCallsFile(FileWriter writer, String fileName, List<cellsGroup> cells, String[] BAMFiles,
									char flag) {
		try {
			writer.write("Gene Name");
			writer.append(',');
			writer.write("Chromosome");
			writer.append(',');
			writer.write("Start");
			writer.append(',');
			writer.write("End");
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++)
					for (int j = cells.get(i).getFirst(); j <= cells.get(i).getLast(); j++) {
						writer.append(',');
						writer.write("cell" + j);
					}
			} else {
				for (int i = 0; i < BAMFiles.length; i++) {
					writer.append(',');
					writer.write(BAMFiles[i]);
				}
			}
			if (flag == 's' || flag == 'l') {
				writer.append(',');
				writer.write("Total");
			}
			writer.append('\n');
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// creates the interval list file
	public void createIntervalList(Gene gene, int startingPosition, int endingPosition,
								   String bamFolder, String fileName,
								   String outputFolder) {
		try {
			String intervalListName = outputFolder + "/one_gene_" + fileName + ".interval_list";

			String command = "java -jar /gpfs0/tals/projects/Analysis/Amit_code/Shai_work/picard.jar ViewSam I="+ bamFolder +"/"+ fileName + " HEADER_ONLY=true";

			command += " > " + intervalListName;
			String[] cmd = { "/bin/sh", "-c", command };
			Process p = Runtime.getRuntime().exec(cmd);
			p.waitFor();
			Runtime.getRuntime().exec("kill $p");

			File file = new File(intervalListName);
			FileWriter fr = new FileWriter(file, true);
			BufferedWriter br = new BufferedWriter(fr);

			br.write(gene.getChrom());
			br.write("\t");
			br.write(String.format("%d", startingPosition));
			br.write("\t");
			br.write(String.format("%d", endingPosition));
			br.write("\t");
			br.write("-");
			br.write("\t");
			br.write(gene.getName());

			if (!gene.getChrom().toLowerCase().startsWith("chr"))
				System.out.println("The chromosome of the gene" + gene.getName() + "does not start with 'chr'");

			br.close();
			fr.close();

		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}

	}

	// write in the junctions file the junctions for every gene
	public void writeJunction(FileWriter writer, Gene gene, List<cellsGroup> cells, int numOfCalls, char flag) {
		try {
			for (int i = 0; i < gene.getJunctions().length; i++) {
				Junction junc = (Junction) gene.getJunctions()[i];
				writer.write(gene.getName());
				writer.append(',');
				writer.write(gene.getChrom());
				writer.append(',');
				if (flag == 's' || flag == 'l') {
					if (junc.getStartCounter() != 0 && junc.getEndCounter() != 0) {
						writer.write(String.format("%d", junc.getStartCounter()));
						writer.write(";");
						writer.write(String.format("%d", junc.getEndCounter()));
					} else if (junc.getStartCounter() != 0)
						writer.write(String.format("%d", (junc.getStartCounter())));
					else
						writer.write(String.format("%d", (junc.getEndCounter())));
					writer.append(',');
				}
				writer.write(String.format("%d", junc.getStart()));
				writer.append(',');
				writer.write(String.format("%d", junc.getEnd()));
				writer.append(',');
				if (flag == 's' || flag == 'l') {
					for (int j = 0; j < cells.size(); j++) {
						for (int a = 0; a < cellsGroup.rangeSize(cells, j); a++) {
							if (numOfCalls == -1)
								writer.write(0);
							else if (junc.getCells().get(j)[a] == 0)
								writer.write("   ");
							else
								writer.write(String.format("%d", junc.getCells().get(j)[a]));
							writer.append(',');
						}
					}

					writer.write(String.format("%d", junc.sumCells()));
				} else {
					for (int j = 0; j < junc.getBamFiles().length; j++) {
						if (junc.getBamFiles()[j] == 0)
							writer.write("   ");
						else
							writer.write(junc.getBamFiles()[j] + "");
						writer.append(',');
					}
				}
				writer.append(System.lineSeparator());

			}

		} catch (IOException e) {
			e.printStackTrace();

		}
	}

	// write for every gene in which cell and exon the calls fall
	public void writeExon(FileWriter writer, Gene gene, List<cellsGroup> cells, int numOfCalls, char flag) {
		try {
			for (int i = 0; i < gene.getExons().length; i++) {
				Exon exon = (Exon) gene.getExons()[i];
				writer.write(gene.getName());
				writer.append(',');
				writer.write(gene.getChrom());
				writer.append(',');
				writer.write(String.format("%d", exon.getStart()));
				writer.append(',');
				writer.write(String.format("%d", exon.getEnd()));
				writer.append(',');
				if (flag == 's' || flag == 'l') {
					for (int j = 0; j < cells.size(); j++) {
						for (int a = 0; a < cellsGroup.rangeSize(cells, j); a++) {
							if (numOfCalls == -1)
								writer.write(0);
							else if (exon.getCells().get(j)[a] == 0)
								writer.write("   ");
							else
								writer.write(String.format("%d", exon.getCells().get(j)[a]));
							writer.append(',');
						}
					}
					writer.write(String.format("%d", exon.sumCells()));
				} else {
					for (int a = 0; a < exon.getBamFiles().length; a++) {
						if (exon.getBamFiles()[a] == 0)
							writer.write("   ");
						else
							writer.write(exon.getBamFiles()[a] + "");
						writer.append(',');
					}

				}
				writer.append(System.lineSeparator());

			}
			// skippingExon(41310655, 41310762, writer2, gene, cells);
			// writePercentage(41310655, 41310762, writer2, gene, cells);
		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	// write in the introns file the introns for every gene
	public void writeIntron(FileWriter writer, Gene gene, List<cellsGroup> cells, int numOfCalls, char flag) {
		try {
			for (int i = 0; i < gene.getIntrons().length; i++) {
				Intron intron = (Intron) gene.getIntrons()[i];
				writer.write(gene.getName());
				writer.append(',');
				writer.write(gene.getChrom());
				writer.append(',');
				writer.write(String.format("%d", intron.getStart()));
				writer.append(',');
				writer.write(String.format("%d", intron.getEnd()));
				writer.append(',');
				if (flag == 's' || flag == 'l') {
					for (int j = 0; j < cells.size(); j++) {
						for (int a = 0; a < cellsGroup.rangeSize(cells, j); a++) {
							if (numOfCalls == -1)
								writer.write(0);
							else if (intron.getCells().get(j)[a] == 0)
								writer.write("   ");
							else
								writer.write(String.format("%d", intron.getCells().get(j)[a]));
							writer.append(',');
						}
					}
					writer.write(String.format("%d", intron.sumCells()));
				} else {
					for (int a = 0; a < intron.getBamFiles().length; a++) {
						if (intron.getBamFiles()[a] == 0)
							writer.write("   ");
						else
							writer.write(intron.getBamFiles()[a] + "");
						writer.append(',');
					}
				}
				writer.append(System.lineSeparator());

			}
		} catch (IOException e) {
			e.printStackTrace();

		}
	}

	// write for every gene in which cell and exon the calls fall
	public void writeUnknownJunction(FileWriter writer, Gene gene, List<cellsGroup> cells, int geneStart, int geneEnd,
									 int numOfCalls, char flag) {
		try {

			writer.write(gene.getName());
			writer.append(',');
			writer.write(gene.getChrom());
			writer.append(',');
			writer.write(String.format("%d", geneStart));
			writer.append(',');
			writer.write(String.format("%d", geneEnd));
			writer.append(',');
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++) {
					for (int a = 0; a < cellsGroup.rangeSize(cells, i); a++) {
						if (numOfCalls == -1)
							writer.write(0);
						else if (gene.getCells().get(i)[a] == 0)
							writer.write("   ");
						else
							writer.write(String.format("%d", gene.getUnknowncells().get(i)[a]));
						writer.append(',');
					}
				}
				writer.write(String.format("%d", gene.sumUnknownJunctionCalls()));
			} else {
				for (int i = 0; i < gene.getUnknownBF().length; i++) {
					if (gene.getUnknownBF()[i] == 0)
						writer.write("   ");
					else
						writer.write(gene.getUnknownBF()[i] + "");
					writer.append(',');
				}
			}
			writer.append(System.lineSeparator());
		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	public void writeGeneCell(FileWriter writer, Gene gene, List<cellsGroup> cells, int geneStart, int geneEnd,
							  int totalCalls, int numOfCalls, char flag) {
		try {
			writer.write(gene.getName());
			writer.append(',');
			writer.write(gene.getChrom());
			writer.append(',');
			writer.write(String.format("%d", geneStart));
			writer.append(',');
			writer.write(String.format("%d", geneEnd));
			writer.append(',');
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++) {
					for (int a = 0; a < cellsGroup.rangeSize(cells, i); a++) {
						if (numOfCalls == -1)
							writer.write(0);
						else if (gene.getCells().get(i)[a] == 0)
							writer.write("   ");
						else
							writer.write(String.format("%d", gene.getCells().get(i)[a]));
						writer.append(',');
					}
				}
				writer.write(String.format("%d", totalCalls));
			} else {
				for (int i = 0; i < gene.getTotalCallsBF().length; i++) {
					if (gene.getTotalCallsBF()[i] == 0)
						writer.write("   ");
					else
						writer.write(gene.getTotalCallsBF()[i] + "");
					writer.append(',');

				}
			}
			writer.append(System.lineSeparator());

		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	public void writedifCigar(FileWriter writer, Gene gene, List<cellsGroup> cells, int geneStart, int geneEnd,
							  int numOfCalls, char flag) {
		try {

			writer.write(gene.getName());
			writer.append(',');
			writer.write(gene.getChrom());
			writer.append(',');
			writer.write(String.format("%d", geneStart));
			writer.append(',');
			writer.write(String.format("%d", geneEnd));
			writer.append(',');
			if (flag == 's' || flag == 'l') {
				for (int i = 0; i < cells.size(); i++) {
					for (int a = 0; a < cellsGroup.rangeSize(cells, i); a++) {
						if (numOfCalls == -1)
							writer.write(0);
						else if (gene.getDifCigarCells().get(i)[a] == 0)
							writer.write("   ");
						else
							writer.write(String.format("%d", gene.getDifCigarCells().get(i)[a]));
						writer.append(',');
					}
				}
				writer.write(String.format("%d", gene.sumDifCigarCalls()));
			} else {
				for (int i = 0; i < gene.getDifCigarBF().length; i++) {
					if (gene.getDifCigarBF()[i] == 0)
						writer.write("   ");
					else
						writer.write(gene.getDifCigarBF()[i] + "");
					writer.append(',');
				}

			}
			writer.append(System.lineSeparator());

		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	private void skippingExon(int exonStart, int exonEnd, FileWriter file, Gene gene, List<cellsGroup> cells) {
		try {
			float counter = 0;
			float numOfCells = 0;
			for (int i = 0; i < 4; i++)
				file.append(',');
			for (int i = 0; i < gene.getExons().length; i++) {
				if (gene.getExons()[i].getStart() == exonStart && gene.getExons()[i].getEnd() == exonEnd) {
					for (int j = 0; j < cells.size(); j++) {
						for (int a = 0; a < cellsGroup.rangeSize(cells, j); a++) {
							numOfCells++;
							if (((Exon) gene.getExons()[i]).getCells().get(j)[a] > 5) {
								file.write("X");
								counter++;
							} else
								file.write("   ");
							file.append(',');
						}
					}
					float t2 = (float) counter / numOfCells;
					file.write(String.format("%f", t2 * 100));
					file.append(System.lineSeparator());
				}
			}
		} catch (IOException e) {
			e.getMessage();
		}

	}

	private void writePercentage(int exonStart, int exonEnd, FileWriter file, Gene gene, List<cellsGroup> cells) {
		try {
			for (int i = 0; i < 4; i++)
				file.append(',');
			for (int i = 0; i < gene.getExons().length; i++) {
				if (gene.getExons()[i].getStart() == exonStart && gene.getExons()[i].getEnd() == exonEnd) {
					float totalCalls = (float) ((Exon) gene.getExons()[i]).sumCells();
					for (int j = 0; j < cells.size(); j++) {
						for (int a = 0; a < cellsGroup.rangeSize(cells, j); a++) {
							float temp = (float) ((Exon) gene.getExons()[i]).getCells().get(j)[a];
							float prec = (float) temp / totalCalls;
							file.write(String.format("%f", (prec * 100)));
							file.append(',');
						}
					}
					file.append(System.lineSeparator());
				}
			}
		} catch (IOException e) {
			e.getMessage();
		}

	}

	// close files
	public void closeFiles(FileWriter[] files) throws IOException {
		try {
			for (FileWriter f : files)
				f.close();
		} catch (IOException e) {
			e.getMessage();
		}
	}

	public FileWriter createRunLaterFile(String outputPathFolder, String runLaterFileName) {
		FileWriter bw = null;
		try {

			bw = new FileWriter(outputPathFolder + runLaterFileName, true);
			File directory = new File(outputPathFolder + "/LG");
			File[] fList = directory.listFiles();
			if (fList != null) {
				for (int i = 0; i < fList.length; i++) {
					if (fList[i].isFile()) {
						BufferedReader br = new BufferedReader(new FileReader(fList[i]));
						String line = "";
						while ((line = br.readLine()) != null)
							bw.write(line + '\n');
						br.close();
						fList[i].delete();
					}
				}
			}
			bw.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return bw;
	}

	//read transcripts file into hashmap
	public HashMap<String, Gene> readTranscriptsFileIntoHashMap(BufferedReader transcriptFile, List<cellsGroup> cells, int bamfilesNumber) {
		HashMap<String, Gene> genes = new HashMap<String, Gene>();
		String line;
		try {
			while ((line = transcriptFile.readLine()) != null) {
				String[] splitLine = line.split("\t");
				if (!genes.containsKey(splitLine[0]) && splitLine.length > 2) {
					if (!splitLine[1].toLowerCase().startsWith("chr"))
						System.out.println("The chromosome of " + splitLine[0] + " does not start with 'chr'");
					Gene gene = new Gene(splitLine[0], splitLine[1], cells, new LinkedList<String>(), new LinkedList<String>(), bamfilesNumber);
					genes.put(splitLine[0], gene);
				}
				genes.get(splitLine[0]).getExonStarts().add(splitLine[2]);
				genes.get(splitLine[0]).getExonEnds().add(splitLine[3]);

			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return genes;

	}

	public void writeGenesHashSetIntoFile(String fileName, HashMap<String, Gene> genes) {
		try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
				for(String geneName : genes.keySet()) {
					genes.get(geneName).writeGeneToFile(bw);
					bw.write('\n');
				}
				bw.close();
			}

		catch (IOException e){
			e.printStackTrace();
		}
	}

}
