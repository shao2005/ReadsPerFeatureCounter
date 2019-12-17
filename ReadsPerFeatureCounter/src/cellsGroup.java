import java.util.ArrayList;
import java.util.List;

public class cellsGroup {
	private int first;
	private int last;

	public int getFirst() {
		return first;
	}

	public cellsGroup(int first, int last) {
		super();
		this.first = first;
		this.last = last;
	}

	public void setFirst(int first) {
		this.first = first;
	}

	public int getLast() {
		return last;
	}

	public void setLast(int last) {
		this.last = last;
	}

	public String toString() {
		return this.getFirst() + "-" + this.getLast();
	}

	// find on which of the cell ranges the call is falling
	public static int findCellsRange(List<cellsGroup> cells, int cell) {
		for (int i = 0; i < cells.size(); i++)
			if (cell >= cells.get(i).getFirst() && cell <= cells.get(i).getLast())
				return i;
		return -1;
	}

	// find the range size
	public static int rangeSize(List<cellsGroup> cells, int i) {
		return cells.get(i).getLast() - cells.get(i).getFirst() + 1;
	}

	// create cell list according to the given ranges
	public static List<int[]> createCellsList(List<cellsGroup> cells) {
		List<int[]> cellsList = new ArrayList<int[]>();
		for (int i = 0; i < cells.size(); i++) {
			cellsList.add(new int[cells.get(i).getLast() + cells.get(i).getFirst() + 1]);
		}
		return cellsList;
	}

}
