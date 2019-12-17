public class Pair {
	
	private String sign;
	private String name;
	
	
	public Pair(String sign, String name) {
		super();
		this.sign = sign;
		this.name = name;
	}
	
	public String getSign() {
		return sign;
	}
	public void setSign(String sign) {
		this.sign = sign;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}

	@Override
	public String toString() {
		return "Pair ["+sign + ",  "+ name+ "]";
	}
	
	

}
