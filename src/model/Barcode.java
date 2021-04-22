package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class Barcode 
{
	public String first;
	public String second;
	public String name;
	public String id;
	public int count = 0;
	
	public static ArrayList<Barcode> readBarcodeFile(File filename)
	{
		ArrayList<Barcode> bc = new ArrayList<Barcode>();
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = br.readLine(); // One barcode per line
			while(line != null)
			{
				String[] tokens = line.split("\t");
				Barcode b = new Barcode();
				b.name = tokens[0].trim();
				b.id = tokens[1].trim();
				b.first = tokens[2].toUpperCase().trim();
				if(b.first.length() != Parameters.lBC1) new ErrorMessage("Barcode " + b.first + " is not of the correct length (" + Parameters.lBC1 + " )");
				b.second = tokens[3].toUpperCase().trim();
				if(b.second.length() != Parameters.lBC2) new ErrorMessage("Barcode " + b.second + " is not of the correct length (" + Parameters.lBC2 + " )");
				bc.add(b);
				line = br.readLine();
			}
			br.close();
		}
		catch(IOException ioe)
		{
			ioe.printStackTrace();
		}
		return bc;
	}
	
	public static HashSet<String> getUniqueBarcodes(ArrayList<Read> list) // TODO Handle sequencing mismatches
	{
		HashSet<String> barcodes = new HashSet<>();
		for(Read r:list)
		{
			barcodes.add(r.barcode);
		}
		return barcodes;
	}
	
	public static String construct(String a, int repeat)
	{
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < repeat; i++) sb.append(a);
		return sb.toString();
	}
}
