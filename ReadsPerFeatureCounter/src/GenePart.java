import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

public abstract class GenePart {
	private String geneName;
	private int start;
	private int end;
	private Hashtable<String, Integer> callsID;
	private List<int[]> cells;
	private int[] bamFiles;

	public GenePart(String geneName, int start, int end, List<int[]> cells, int numOfBamFiles) {
		super();
		this.geneName = geneName;
		this.start = start;
		this.end = end;
		callsID = new Hashtable<>();
		this.cells = cells;
		bamFiles = new int[numOfBamFiles];
	}

	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public  Hashtable<String, Integer> getCallsID() {
		return callsID;
	}

	public void setCallsID( Hashtable<String, Integer> callsID) {
		this.callsID = callsID;
	}

	public List<int[]> getCells() {
		return cells;
	}

	public void setCells(List<int[]> cells) {
		this.cells = cells;
	}

	public int[] getBamFiles() {
		return bamFiles;
	}

	public void setBamFiles(int[] bamFiles) {
		bamFiles = bamFiles;
	}

	public int sumCells() {
		int sum = 0;
		for (int i = 0; i < cells.size(); i++)
			for (int j = 0; j < cells.get(i).length; j++)
				sum += cells.get(i)[j];
		return sum;
	}

	@Override
	public String toString() {
		return "[" + start + "-" + end + "]";
	}

}
