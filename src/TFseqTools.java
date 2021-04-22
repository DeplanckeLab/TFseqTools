import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import model.Barcode;
import model.ErrorMessage;
import model.Parameters;
import model.Read;
import model.UMI;
import tools.Logger;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class TFseqTools 
{	
	public static void main(String[] args)
	{
		System.out.println("TFseqTools v"+Parameters.currentVersion);
		if(args.length < 1) Parameters.printHelp();
		else
		{
			// Parsing args
			String[] argsParsed = new String[args.length - 1];
			for(int i = 0; i < args.length - 1; i++) argsParsed[i] = args[i + 1];
			
			switch(args[0])
			{
				case "Counter":
					Parameters.loadCounter(argsParsed);
					Logger.init(Parameters.logFile);
					
					System.out.println("\n-1- |  Reading TF barcodes");
					Parameters.bc = Barcode.readBarcodeFile(Parameters.inputTFFile);
					System.out.println(Parameters.bc.size() + " barcodes found in barcode file.");
					
					System.out.println("\n-2- |  Reading BAM file");
					HashMap<String, Barcode> mappedReads = AnalyzeAlignedBAM.readR2BAM(Parameters.inputBAMFileR2);
					System.out.println(Parameters.nbReads + " total reads in BAM file.");
					System.out.println(Parameters.unmapped + " unmapped reads in BAM file.");
					System.out.println(Parameters.notUnique + " not unique alignments.");
					System.out.println(Parameters.tooLowAQUAL + " too Low AQual.");
					System.out.println(Parameters.tooLowSQUAL + " too Low SQual.");
					System.out.println(mappedReads.size() + " reads are consistently mapping to existing TFs.");
					System.out.println((Parameters.nbReads - mappedReads.size()) + " reads are not mapping to TFs\t(" + Parameters.myFormatter.format(((Parameters.nbReads - mappedReads.size()) / (float)Parameters.nbReads) * 100) + "%)");
					System.out.println(Parameters.overlapBC1 + " reads overlap with barcode 1 only\t(" + Parameters.myFormatter.format((Parameters.overlapBC1 / (float)Parameters.nbReads) * 100) + "%)");
					System.out.println(Parameters.overlapBC2 + " reads overlap with barcode 2 only\t(" + Parameters.myFormatter.format((Parameters.overlapBC2 / (float)Parameters.nbReads) * 100) + "%)");
					System.out.println(Parameters.overlapBoth + " reads overlap with both barcodes\t(" + Parameters.myFormatter.format((Parameters.overlapBoth / (float)Parameters.nbReads) * 100) + "%)");
					
					System.out.println("\n-3- |  Reading R1 fastq file");
					ArrayList<Read> finalReads = Utils.readR1Fastq(mappedReads);
					System.out.println(Parameters.nbReads + " total reads in FASTQ file.");
					System.out.println(finalReads.size() + " R2 TF reads were matching R1 fastq file");
					HashSet<String> uniqueBarcodes = Barcode.getUniqueBarcodes(finalReads);
					System.out.println(uniqueBarcodes.size() + " unique CELL barcodes were found.");
					
					Logger.close();
					
					// Build indexes for filling the count/UMI matrix
					HashMap<String, Integer> barcodeIndexes = new HashMap<>();
					int index = 0;
					for(String bc:uniqueBarcodes) {barcodeIndexes.put(bc, index); index++;}
					HashMap<String, Integer> TFIndexes = new HashMap<>();
					index = 0;
					for(Barcode bc:Parameters.bc) {TFIndexes.put(bc.id, index); index++;}
					
					// Fill the count matrix
					int[][] countMatrix = new int[Parameters.bc.size()][uniqueBarcodes.size()];
					for(Read r:finalReads) countMatrix[TFIndexes.get(r.tfBarcode.id)][barcodeIndexes.get(r.barcode)]++;
					
					// Fill the UMI matrix
					UMI[][] umiMatrix = new UMI[Parameters.bc.size()][uniqueBarcodes.size()];
					for(int i = 0; i < umiMatrix.length; i++) for(int j = 0; j < umiMatrix[i].length; j++) umiMatrix[i][j] = new UMI(); 
					for(Read r:finalReads) umiMatrix[TFIndexes.get(r.tfBarcode.id)][barcodeIndexes.get(r.barcode)].addUMI(r.UMI);
					
					try
					{
						BufferedWriter results = new BufferedWriter(new FileWriter(Parameters.outputFolder + "Results.Matrix.txt"));
						results.write("TFName\tTFId");
						for(String bc:uniqueBarcodes) results.write("\t" + bc);
						results.write("\n");
						for(Barcode b:Parameters.bc) 
						{
							results.write(b.name + "\t" + b.id);
							for(String bc:uniqueBarcodes)
							{
								results.write("\t" + countMatrix[TFIndexes.get(b.id)][barcodeIndexes.get(bc)]);
							}
							results.write("\n");
						}
						results.close();
					}
					catch(IOException ioe)
					{
						new ErrorMessage(ioe.getMessage());
					}
					
					try
					{
						BufferedWriter results = new BufferedWriter(new FileWriter(Parameters.outputFolder + "Results.Matrix.UMI.txt"));
						results.write("TFName\tTFId");
						for(String bc:uniqueBarcodes) results.write("\t" + bc);
						results.write("\n");
						for(Barcode b:Parameters.bc) 
						{
							results.write(b.name + "\t" + b.id);
							for(String bc:uniqueBarcodes)
							{
								results.write("\t" + umiMatrix[TFIndexes.get(b.id)][barcodeIndexes.get(bc)].getCorrectedSize());
							}
							results.write("\n");
						}
						results.close();
					}
					catch(IOException ioe)
					{
						new ErrorMessage(ioe.getMessage());
					}
					
					break;
				default:
					Parameters.printHelp();
					new ErrorMessage("The tool '"+ args[0] +"' is not implemented. Please use one of the following: [Counter].");
			}
		}
	}
	

}

