package tools;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import model.Barcode;
import model.ErrorMessage;
import model.Parameters;
import model.Read;

public class Utils 
{
	private static Random rand = new Random();
	
	/**
	 * Read FastQ file containing the mapping between read name, barcode and UMI.
	 * @param fastQ input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static BufferedReader readFastq(File fastQ)
	{
		if(fastQ.getAbsolutePath().endsWith(".fastq") || fastQ.getAbsolutePath().endsWith(".fq"))
		{
			try
			{
				return new BufferedReader(new FileReader(fastQ));
			}
			catch(FileNotFoundException fnfe)
			{
				new ErrorMessage(fnfe.getMessage());
			}
		}
		else if(fastQ.getAbsolutePath().endsWith(".fastq.gz") || fastQ.getAbsolutePath().endsWith(".fq.gz"))
		{
			try
			{
				return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastQ))));
			}
			catch(IOException ioe)
			{
				new ErrorMessage(ioe.getMessage());
			}
		}
		else
		{
			new ErrorMessage("The extension of the FastQ file is not recognized : " + fastQ.getAbsolutePath() + "\nIt should be '.fastq', '.fq', '.fq.gz' or '.fastq.gz'");
		}
		return null;
	}
	
	/**
	 *  Reading reads barcodes/UMI from the R1 fastq file to map the UMI/barcode with the read name (lost after alignment)
	 * @throws Exception Yes I know...
	 */
	public static ArrayList<Read> readR1Fastq(HashMap<String, Barcode> mappedReads)
	{
		ArrayList<Read> finalReads = new ArrayList<Read>();
		
		System.out.println("\nReading reads barcodes/UMI from the R1 fastq file...");
		BufferedReader br = Utils.readFastq(Parameters.inputFastQFileR1);
		Long start = System.currentTimeMillis();

		Parameters.nbReads = 0;
		
		Read read = Utils.nextRead(br);
		while(read != null)
		{
			Parameters.nbReads++;
			if(Parameters.nbReads %Parameters.chunkSize == 0) System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			Barcode bc = mappedReads.get(read.name);
			if(bc != null)
			{
				read.tfBarcode = bc;
				finalReads.add(read);
			}
			read = Utils.nextRead(br);
		}
		
		Utils.close(br);
				
		System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");		
		
		return finalReads;
	}
	
	/**
	 * Read FastQ file containing the mapping between read name, barcode, UMI
	 * @param fastQ input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static Read nextRead(BufferedReader br)
	{
		Read read = new Read();
		try
		{
			String header = br.readLine(); // Get first line = READNAME + INDEX
			if(header == null) return null;
			if(!header.startsWith("@")) new ErrorMessage("R1 fastq file has formatting issues");
			read.name = header.substring(1, header.indexOf(" "));
			String indexes = br.readLine().trim(); // Get second line = READ
			if(indexes.length() != Parameters.lengthBarcode) // Checking the length is OK
			{
				br.close();
				System.err.println("Error while parsing R1 FastQ file");
				System.err.println("Read found in FastQ has length " + indexes.length() + " while barcode pattern has length " + Parameters.lengthBarcode);
				System.exit(-1);
			}
			try
			{
				if(Parameters.l1 != -1) read.barcode = indexes.substring(Parameters.barcodeRange[0], Parameters.barcodeRange[1]); // If there is a barcode to look for
				if(Parameters.UMILength != -1) read.UMI = indexes.substring(Parameters.UMIRange[0], Parameters.UMIRange[1]);
			}
			catch(IndexOutOfBoundsException boe)
			{
				br.close();
				System.err.println("Error while parsing R1 FastQ file");
				System.err.println("Index found in FastQ "+ indexes +" is not according to the pattern: " + Parameters.barcodePattern);
				System.exit(-1);
			}
			br.readLine(); // Third line = we don't care
			br.readLine();  // Fourth line = qualities = we don't care
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		return read;
	}
	
	/**
	 * Parsing the barcode pattern and check validity
	 * @throws Exception Yes I know...
	 */
	public static void patterning()
	{
		Parameters.lengthBarcode = 0;
		for(int j = 0; j < Parameters.barcodePattern.length(); j++) 
		{
			char c = Parameters.barcodePattern.charAt(j);
			switch(c)
			{
				case 'B':
					Parameters.barcodeRange[0] = Parameters.lengthBarcode;
					Parameters.lengthBarcode += Parameters.l1;
					Parameters.barcodeRange[1] = Parameters.lengthBarcode;
					break;					
				case 'U':
					Parameters.UMIRange[0] = Parameters.lengthBarcode;
					Parameters.lengthBarcode += Parameters.UMILength;
					Parameters.UMIRange[1] = Parameters.lengthBarcode;
					break;
				case '?':
					Parameters.lengthBarcode++;
					break;
				default:
					new ErrorMessage(c+" does not correspond to any authorized pattern character (?, B, U). Aborted.");
			}
		}
		System.out.println("\tAccording to barcode pattern, reads of R1 FastQ file should contain "+Parameters.lengthBarcode+" characters.");
	}
	
	public static float min(float a, float b, float c) 
	{
		return Math.min(Math.min(a, b), c);
	}
	
	public static void setSeed(int seed)
	{
		rand = new Random(seed);
	}
	
	public static void listdirs(String directoryName, ArrayList<File> files) 
	{
	    File directory = new File(directoryName);
	    File[] fList = directory.listFiles();
	    for (File file : fList) if (file.isDirectory()) files.add(file);
	}
	
	public static void listfiles(String directoryName, ArrayList<File> files) 
	{
	    File directory = new File(directoryName);
	    File[] fList = directory.listFiles();
	    for (File file : fList) if (file.isFile()) files.add(file);
	}
    
	public static Random getRandomGenerator()
	{
		return rand;
	}
	
	public static double mean(double[] data)
	{
		double sum = 0.0;
		for(double a : data) sum += a;
		return sum/data.length;
	}
	
	public static double mean(ArrayList<Integer> data)
	{
		double sum = 0.0;
		for(Integer a : data) sum += a;
		return sum/(double)data.size();
	}
	
	public static int sum(int[] array)
	{
		int sum = 0;
		for(int i:array) sum+=i;
		return sum;
	}
	
	public static double mean(int[] data)
	{
		double sum = 0.0;
		for(double a : data) sum += a;
		return sum/data.length;
	}

	public static double median(int[] data)
	{
		Arrays.sort(data);
		if (data.length % 2 == 0) return ((double)data[data.length/2] + (double)data[data.length/2 - 1])/2;
		return (double) data[data.length/2];
	}
	
	public static double var(double[] data, double mean)
	{
		double var = 0;
		for(double a :data) var += (a-mean)*(a-mean);
		return var/(data.length-1);
	}
	
	public static double var(int[] data, double mean)
	{
		double var = 0;
		for(double a :data) var += (a-mean)*(a-mean);
		return var/(data.length-1);
	}
	
	public static double sd(double[] data, double mean)
	{
		double var = var(data, mean);
		return Math.sqrt(var);
	}
	
	public static double sd(int[] data, double mean)
	{
		double var = var(data, mean);
		return Math.sqrt(var);
	}
	
	public static double var(double[] data)
	{
	    double mean = 0;
	    double M2 = 0;
	    if(data.length < 2) return 0;
	    for (int i = 0; i < data.length; i++) 
	    {
	        double delta = data[i] - mean;
	        mean = mean + delta/(i+1);
	        M2 = M2 + delta*(data[i] - mean);
		}
	    return M2/(data.length - 1);
	}
	
	public static double sd(double[] data)
	{
		double var = var(data);
		return Math.sqrt(var);
	}
	
	public static double cv(double[] data) // Coefficient of Variation
	{
		double mu = mean(data);
		double theta = sd(data);
		return theta / mu;
	}
	
	public static double[] toArray(ArrayList<Double> data)
	{
		double[] res = new double[data.size()];
		for (int i = 0; i < res.length; i++) res[i] = data.get(i);
		return res;
	}
	
	public static void close(BufferedReader br)
	{
		if(br == null) new ErrorMessage("This file handle is not initialized.");
		try
		{
			br.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
	}
	
    public static double cov(double[] x, double[] y)
    {
        double result = 0d;
        int length = x.length;
        double xMean = mean(x);
        double yMean = mean(y);
        for (int i = 0; i < length; i++) 
        {
        	double xDev = x[i] - xMean;
        	double yDev = y[i] - yMean;
        	result += (xDev * yDev - result) / (i + 1);
        }
        return result * ((double) length / (double)(length - 1)); // Bias correction
    }
    
    public static double cov(double[] x, double[] y, double xMean, double yMean)
    {
        double result = 0d;
        int length = x.length;
        for (int i = 0; i < length; i++) 
        {
        	double xDev = x[i] - xMean;
        	double yDev = y[i] - yMean;
        	result += (xDev * yDev - result) / (i + 1);
        }
        return result * ((double) length / (double)(length - 1)); // Bias correction
    }
	
	public static String[] sortD(Map<String, Double> map)
	{
		List<Double> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static boolean contains(String[] array, String value)
	{
		for(String val:array) if(val.equals(value)) return true;
		return false;
	}
	
	public static int hammingDistance(String a, String b) 
	{
		// Check if we can compute the distance
		if (a == null || b == null) new ErrorMessage("Strings must not be null");
		int l = a.length();
		if (l != b.length()) new ErrorMessage("Strings must have the same length");
		
		// Compute the distance
		int distance = 0;
		for(int i = 0; i < l; i++) if (a.charAt(i) != b.charAt(i)) distance++;
		
		return distance;
	}

	public static String[] sortI(Map<String, Integer> map)
	{
		List<Integer> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static String[] sortL(Map<String, Long> map)
	{
		List<Long> values = new ArrayList<>(map.values());
		Collections.sort(values, Collections.reverseOrder());
		String[] sortedIndexes = new String[values.size()];
		for(String key:map.keySet()) 
		{
			int index = values.indexOf(map.get(key));
			sortedIndexes[index] = key;
			values.set(index, null);
		}
		return sortedIndexes;
	}
	
	public static String[] sortKeys(Set<String> keySet)
	{
		String[] keys = keySet.toArray(new String[keySet.size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String toReadableTime(long ms)
	{
		if(ms < 1000) return ""+ms+" ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if(s < 60) return s+" s "+ms+" ms";
		long mn = s / 60;
		s = s % 60;
		if(mn < 60) return mn +" mn "+s+" s "+ms+" ms";
		long h = mn / 60;
		mn = mn % 60;
		if(h < 24) return h +" h "+ mn +" mn "+s+" s "+ms+" ms";
		long d = h / 24;
		h = h % 24;
		return d+ " d " + h +" h "+ mn +" mn "+s+" s "+ms+" ms";
	}
}
