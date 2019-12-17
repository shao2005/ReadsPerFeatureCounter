import java.util.*;

public class Junction extends GenePart {
	private int length;
	private int startCounter;
	private int endCounter;

	public Junction(String geneName, int start, int end, List<int[]> cells, int length, int startCounter,
			int endCounter, int numOfBamFiles) {
		super(geneName, start, end, cells, numOfBamFiles);
		this.length = length;
		this.startCounter = startCounter;
		this.endCounter = endCounter;
	}

	public Junction(Junction junction) {
		super(junction.getGeneName(), junction.getStart(), junction.getEnd(), junction.getCells(),
				junction.getBamFiles().length);
		this.length = junction.length;
		this.startCounter = junction.startCounter;
		this.endCounter = junction.endCounter;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public int getStartCounter() {
		return startCounter;
	}

	public void setStartCounter(int startCounter) {
		this.startCounter = startCounter;
	}

	public int getEndCounter() {
		return endCounter;
	}

	public void setEndCounter(int endCounter) {
		this.endCounter = endCounter;
	}

}
