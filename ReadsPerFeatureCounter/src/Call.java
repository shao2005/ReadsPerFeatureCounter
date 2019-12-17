public class Call {
	private String[] call;
	private String ID;
	private int start;
	private String cigar;
	private int cell;
	
	public Call(String[] call, String ID,int start, String cigar, int cell) {
		super();
		this.ID=ID;
		this.call=call;
		this.start = start;
		this.cigar = cigar;
		this.cell = cell;		
	}
	
	public String[] getCall() {
		return call;
	}

	public void setCall(String[] call) {
		this.call = call;
	}

	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
	}

	public Call(String[] call){
		try {
			this.call = call;
			this.ID = call[0];
			this.start = Util.convertToInteger(call[3]);
			this.cigar = call[5];
			int cellColumn = findCell();
			this.cell = Util.convertToInteger(call[cellColumn].substring(5));
		}
		catch(Exception e){
			e.printStackTrace();
		}

	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public String getCigar() {
		return cigar;
	}

	public void setCigar(String cigar) {
		this.cigar = cigar;
	}

	public int getCell() {
		return cell;
	}

	public void setCell(int cell) {
		this.cell = cell;
	}
	
	//find the place of the cell in the call
	public int findCell(){
		int cellColumn=0;
		for(int a=0; a<call.length && cellColumn==0; a++){
			if(call[a].contains("xi:"))
				cellColumn=a;
		}
		return cellColumn;
	}
	
	public String toString (){
		return ID+'\t'+start+'\t'+cigar+'\t'+cell;
	}
	
	

}
