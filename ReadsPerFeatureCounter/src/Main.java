import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
 
public class Main {

	static Processing obj = new Processing();

	public static void main(String[] args) throws IOException {
		if(args[0].equals("sc"))
			singleCellRun(args);
		else if(args[0].equals("pop"))
			populationRun(args);
		else if(args[0].equals("trans"))
			createMouseTranscriptsFile(args);

	}

	public static void createMouseTranscriptsFile(String[] args) {
		String workingDirectory = args[1];
		String gtfFile = args[2];
		String transcriptsFileName = args[3];
//		String transcriptsFileName="mouseGenesTest.txt";
		// String gtfFile = "gpfs0/tals/projects/data/Transcriptomes/mm10/mm10_combined_annotated_oneName.gtf";
//		String gtfFile = "Z:\\ess_projects\\data\\Transcriptomes\\mm10\\mm10_combined_annotated_oneName.gtf";
		Hashtable<String, Gene> genes = obj.createGenesListMouse(gtfFile, workingDirectory);
		obj.createGenesFile(transcriptsFileName, genes, workingDirectory);
		System.out.println("finish");
	}

	public static void singleCellRun(String[] args) {
		List<cellsGroup> cells = new ArrayList<>();
		int i = 1;
		String picardPath = args[i++];
		String fileName = args[i++];
		String folderName = args[i++];
		String transFile = args[i++];
		String outputFolder = args[i++];
		String geneToStart = args[i++];
		String uniqueId = args[i++];

		while (i < args.length)
			cells.add(new cellsGroup(Integer.parseInt(args[i++]), Integer.parseInt(args[i++])));

		obj.process(cells, fileName, folderName, transFile, outputFolder, 's', new String[0], null, geneToStart, uniqueId, picardPath);
		System.out.println("finish-partA");
		obj.process(cells, fileName, folderName, transFile, outputFolder + "/LG", 'l', new String[0], null, geneToStart, uniqueId, picardPath);
		System.out.println("finish");
	}

	public static void populationRun(String args[]){
		if(args.length < 7)
			System.out.println("Some parameters are missing");
		String[] BamFiles = new String[args.length - 7];
		int i = 1;
		String picardPath = args[i++];
		String folderName = args[i++];
		String transFile = args[i++];
		String outputFolder = args[i++];
		String geneToStart = args[i++];
		String uniqueId = args[i++];

		int j;
		for (j = 0; j < BamFiles.length; j++)
			BamFiles[j] = args[i + j];

		obj.process(null,BamFiles[0], "", transFile, outputFolder, 'p', BamFiles, folderName, geneToStart, uniqueId, picardPath);
		System.out.println("finish-partA");
		obj.process(null, BamFiles[0], "", transFile, outputFolder + "/LG", 'L', BamFiles, folderName, geneToStart, uniqueId, picardPath);
		System.out.println("finish");
	}


	public void createHumanTranscriptsFile(String workingDirectory, String gtfFile, String transcriptsFileName){
//		String gtfFile = "/gpfs0/tals/projects/data/Transcriptomes/human_hg38/ensembl_ucsc_refseq.sorted.gtf";
//		String transcriptsFileName="humanFile_hg38.txt";
		Hashtable<String,Gene> genes = obj.createGenesListHuman(workingDirectory, gtfFile);
		obj.createGenesFile(transcriptsFileName, genes, workingDirectory);
		System.out.println("finish");
	}
}

