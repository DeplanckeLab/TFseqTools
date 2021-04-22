package model;

public class Read implements Comparable<Read>
{
	public String name;
	public String barcode;
	public String UMI;
	public Barcode tfBarcode;
	
	@Override
	public int compareTo(Read r2) 
	{
		return this.name.compareTo(r2.name);
	}
}
