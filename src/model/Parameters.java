package model;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;

import tools.Utils;

enum Strand{NO, YES, REVERSE};

public class Parameters 
{
	public static final String currentVersion = "1.0";
	public static DecimalFormat myFormatter = new DecimalFormat("##.##");
	
	// Barcodes
	public static int startBC1 = 6419; // Pos in fasta genome to align to
	public static int endBC1 = 6429; // Pos in fasta genome to align to
	public static int lBC1 = -1;
	public static int startBC2 = 6457; // Pos in fasta genome to align to
	public static int endBC2 = 6464; // Pos in fasta genome to align to
	public static int lBC2 = -1;
	public static int overlapBC1 = 0;
	public static int overlapBC2 = 0;
	public static int overlapBoth = 0;
	public static ArrayList<Barcode> bc;
	
	// Input parameters
	public static String outputFolder = null;
	public static File logFile = null;
	public static File inputTFFile = null;
	public static File inputFastQFileR1 = null;
	public static File inputBAMFileR2 = null;
	public static long chunkSize = 1000000; // For printing
	
	// Cell barcode in R1
	public static String barcodePattern = "BU";
	public static int hammingDistanceUMI = 0;
	public static int UMILength = -1;
	public static int lengthBarcode = 0;
	public static int l1 = -1;
	public static int[] barcodeRange = new int[2];
	public static int[] UMIRange = new int[2];
	
	// Computed variables
	public static int noFeature = 0;
	public static int notUnique = 0;
	public static int ambiguous = 0;
	public static int mapped = 0;
	public static int unmapped = 0;
	public static int tooLowAQUAL = 0;
	public static int tooLowSQUAL = 0;
	public static long nbReads = 0;
	
	public static void loadCounter(String[] args)
	{
		if(args.length == 0)
		{
			printHelpCounter();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) new ErrorMessage(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "--UMI":
						i++;
						try
						{
							UMILength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorMessage("The '--UMI' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "--tf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputTFFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--tf' option should be followed by the file path containing TF barcodes. " + e.getMessage() + ".");
						}
						break;
					case "--BC":
						i++;
						try
						{
							l1 = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorMessage("The '--BC' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "--nu":
						i++;
						try
						{
							hammingDistanceUMI = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							new ErrorMessage("The '--nu' option should be followed by an Integer. You entered " + args[i]);
						}
						break;
					case "--r1":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputFastQFileR1 = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--r1' option should be followed by R1 FastQ file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--r2":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) new ErrorMessage("No file at path " + args[i]);
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							inputBAMFileR2 = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--r2' option should be followed by aligned BAM file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "--log":
						i++;
						try
						{
							File c = new File(args[i]);
							if(c.exists()) new ErrorMessage("Existing file at path " + args[i]);
							else c.createNewFile();
							if(!c.isFile()) new ErrorMessage(args[i] + " is not a file");
							logFile = c;
						}
						catch(Exception e)
						{
							new ErrorMessage("The '--log' option should be followed by a file path. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					case "-p":
						i++;
						try
						{
							barcodePattern = args[i];
						}
						catch(Exception e)
						{
							new ErrorMessage("The '-p' option should be followed by a valid barcode pattern. " + e.getMessage() + ". You entered " + args[i]);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(inputFastQFileR1 == null)
		{
			new ErrorMessage("Please use '--r1' option to specify R1 FastQ file with barcode/UMI information");
		}
		if(inputBAMFileR2 == null)
		{
			new ErrorMessage("Please use '--r2' option to specify the path of the aligned BAM file");
		}
		if(inputTFFile == null) 
		{
			new ErrorMessage("Please use '--tf' option to specify the path of the TF barcodes to use");
		}
		System.out.println("\n-- Input Parameters --");
		System.out.println("Cell barcodes (in R1 fastq file):");
		System.out.println(" \tPattern = '" + barcodePattern + "'");
		if(barcodePattern.contains("B"))
		{
			if(l1 == -1)
			{
				l1 = 16;
				System.out.println("\tCell barcode sequence length = '"+l1+"' [default to 10x protocol v3, use '--BC' option to change]");
			}
			else System.out.println("\tCell barcode sequence length = '"+l1+"'");
		}
		else 
		{
			if(l1 != -1) new ErrorMessage("Your barcode pattern does not contain any 'B' but you used the '--BC' option to specify a barcode length. Remove the '--BC' option, or change your barcode pattern.");
		}
		if(barcodePattern.contains("U"))
		{
			if(UMILength == -1)
			{
				UMILength = 12;
				System.out.println("\tUMI sequence length = '"+UMILength+"' [default to 10x protocol v3, use '--UMI' option to change]");
			}
			else System.out.println("\tUMI sequence length = '"+UMILength+"'");

		}
		if(UMILength != -1 && !barcodePattern.contains("U"))
		{
			new ErrorMessage("You specified a UMI length but your barcode pattern does not contain 'U', you should specify where to find the UMI in R1");
		}
		Utils.patterning();
		System.out.println("TF barcodes (searched for in R2 bam file):");
		System.out.println("\tReference TF file = " + inputTFFile);
		// Handle barcodes
		Parameters.lBC1 = Parameters.endBC1 - Parameters.startBC1 + 1;
		Parameters.lBC2 = Parameters.endBC2 - Parameters.startBC2 + 1;
		System.out.println("\tBarcode 1 is searched for at pos [" + Parameters.startBC1 + ", " + Parameters.endBC1 + "],  l = " + Parameters.lBC1);
		System.out.println("\tBarcode 2 is searched for at pos [" + Parameters.startBC2 + ", " + Parameters.endBC2 + "],  l = " + Parameters.lBC2);
		
		if(outputFolder == null)
		{
			String path = inputBAMFileR2.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			outputFolder = path;
		}
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
		if(logFile == null) System.out.println("Log File = NONE (specify a log file using option --log)");
		else System.out.println("Log File = " + logFile.getAbsolutePath());
	
		System.out.println("Output folder = " + Parameters.outputFolder + ". Use '-o' option to change.");
	}
	
	private static void printHelpCounter()
	{
		System.out.println("\n-- 'Counter' options --");
		System.out.println("\t--r1 %s \t[Required] Path of R1 FastQ file.");
		System.out.println("\t--r2 %s \t[Required] Path of R2 aligned BAM file [do not need to be sorted or indexed].");
		System.out.println("\t--tf %s \t[Required] File containing known TF barcodes");
		System.out.println("\t-o %s \t\tOutput folder [default = folder of BAM file]");
		System.out.println("\n-- Additional options --");
		System.out.println("\t--log %i \tDetailed log file [default: None]");
		System.out.println("\t--nu %i \tNumber of allowed difference (hamming distance) for two UMIs to be counted only once [default = 0].");
		System.out.println("\t-p %s \t\tCell barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file [default = 'BU', i.e. barcode followed by the UMI].\n\t\t\t\t'B' [Required] is used for specifying the barcode position.\n\t\t\t\t'U' can be used for specifying a UMI value position.\n\t\t\t\t'?' can be used to ignore specific nucleotides.");
		System.out.println("\t--UMI %i \tIf your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI [e.g. 10x run is 10]");
		System.out.println("\t--BC %i \tIf your barcode pattern contains Barcode ('B'), you should specify this parameter as the length of the barcode [e.g. 10x run is 16]");
	}
	
	public static void printHelp()
	{
		System.out.println("\n-- Options --");
		System.out.println("\tCounter\t\tCount occurence of each TF in Drop-seq/10x R2 aligned BAM and R1 fastq");
	}
}
