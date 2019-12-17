import java.util.List;

public class Util {

	//convert string to int
	public static int convertToInteger (String str){
		int ans=0;
		for(int i=0; i<str.length(); i++){		
			char temp=str.charAt(i);
			int digit=temp-'0';
		if(digit>9 || digit<0){
			ans=-1;
			break;
		}
		ans=ans*10+(temp-'0');
			}
			return ans;

	}
	
	//find the minimal start
	public static int findMin(List<List<String>> list)
	{
		int min=Integer.MAX_VALUE;
		for(int i=0; i<list.size(); i++){
			int temp=convertToInteger(list.get(i).get(0));
			if(temp<min)
				min=temp;					
		}
		return min;
	}

	//find the maximal end
	public static int findMax(List<List<String>> list)
	{
		int max=0;
		for(int i=0; i<list.size(); i++){
			int temp=convertToInteger(list.get(i).get(list.get(i).size()-1));
			if(temp>max)
				max=temp;					
		}
		return max;
	}
	
	//find the size of the data
	public static int findSize(String[][] arr){
		int count=0;
		for(int i=0; i<arr.length && arr[i][0]!=null; i++)
			count++;
		return count;
	}
	
	
	

}
