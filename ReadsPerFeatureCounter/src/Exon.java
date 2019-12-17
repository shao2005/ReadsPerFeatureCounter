import java.util.ArrayList;
import java.util.List;

public class Exon extends GenePart {
	private int transcriptionCounter;

	public Exon(String geneName, int start, int end, List<int[]> cells, int transcriptionCounter, int numOfBamFiles) {
		super(geneName, start, end, cells, numOfBamFiles);
		this.transcriptionCounter = transcriptionCounter;
	}

	public int getTranscriptionCounter() {
		return transcriptionCounter;
	}

	public void setTranscriptionCounter(int transcriptionCounter) {
		this.transcriptionCounter = transcriptionCounter;
	}

}
