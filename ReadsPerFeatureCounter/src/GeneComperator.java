import java.util.Comparator;

public class GeneComperator implements Comparator<GenePart>{
	
	 @Override
	    public int compare(GenePart x1, GenePart x2) {
	        return Integer.valueOf(x1.getStart()).compareTo(Integer.valueOf(x2.getStart())) ;     
	    }

}
