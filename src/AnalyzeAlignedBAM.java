import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import model.Barcode;
import model.ErrorMessage;
import model.Parameters;
import model.PolySeq;
import tools.Levenshtein;
import tools.Logger;
import tools.Utils;

public class AnalyzeAlignedBAM 
{	
	/**
	 * Using Picard to read the reads from the BAM file created by the alignment tool
	 * @throws Exception Yes I know...
	 */
	public static HashMap<String, Barcode> readR2BAM(File inputBAMFile)
	{
		HashMap<String, Barcode> result = new HashMap<String, Barcode>();
		Long start = System.currentTimeMillis();
		Logger.write("\nReading the reads from the BAM file...\n");
		try
		{
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = samReaderFactory.open(inputBAMFile);
			SAMRecordIterator it = samReader.iterator();
	
			// Start reading the BAM file
			while(it.hasNext())
			{
				Parameters.nbReads++;
				SAMRecord samRecord = it.next();
				float sequencing_phred = 0;
				for(byte b:samRecord.getBaseQualities()) sequencing_phred += (int)b;
				sequencing_phred /= samRecord.getBaseQualities().length;
				if(samRecord.getSupplementaryAlignmentFlag()) Parameters.notUnique++;
				else if(samRecord.getReadUnmappedFlag()) Parameters.unmapped++;
				else if(samRecord.getMappingQuality() < 10) Parameters.tooLowAQUAL++;
				else if(sequencing_phred < 10) Parameters.tooLowSQUAL++;
				else 
				{
					Barcode bc = overlapsBC(samRecord);
					if(bc != null) result.put(samRecord.getReadName(), bc);
				}
				if(Parameters.nbReads%Parameters.chunkSize == 0) System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			samReader.close();
		}
		catch(IOException ioe)
		{
			new ErrorMessage(ioe.getMessage());
		}
		catch(SAMFormatException sfe)
		{
			new ErrorMessage(sfe.getMessage());
		}
		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		return result;
	}
	
	public static Barcode overlapsBC(SAMRecord samRecord)
	{
		String BC1 = getAlignedStringAtPos(samRecord, Parameters.startBC1, Parameters.endBC1); // position of barcode 1
		String BC2 = getAlignedStringAtPos(samRecord, Parameters.startBC2, Parameters.endBC2); // position of barcode 2
		if(BC1 != null && BC1.equals(Barcode.construct("-", Parameters.lBC1))) BC1 = null; // Construct of length 11
		if(BC2 != null && BC2.equals(Barcode.construct("-", Parameters.lBC2))) BC2 = null; // Construct of length 8
		if(BC1 != null && BC1.equals(Barcode.construct("N", Parameters.lBC1))) BC1 = null; // Construct of length 11
		if(BC2 != null && BC2.equals(Barcode.construct("N", Parameters.lBC2))) BC2 = null; // Construct of length 8
		if(BC1 != null && BC2 != null)
		{
			Levenshtein metric = new Levenshtein();
			Barcode bestMatch1 = null;
			int whichReturn = -1;
			boolean multiple1 = false;
			float max1 = 0;
			String same1 = "Barcodes with same similarity value:";
			for(Barcode b:Parameters.bc)
			{
				float result = metric.compare(BC1, b.first);
				if(result > 0.9) // 0.909... for 1 mismatch 11 bp. The only accepted mismatch rate
				{
					if(result > max1) 
					{
						multiple1 = false;
						same1 = "Barcodes with same similarity value:";
						bestMatch1 = b;
						max1 = result;
					} 
					else if(result == max1) 
					{
						multiple1 = true;
						same1 += "\t" + b.first;
					}
				}
			}

			Barcode bestMatch2 = null;
			float max2 = 0;
			boolean multiple2 = false;
			String same2 = "Barcodes with same similarity value:";
			for(Barcode b:Parameters.bc)
			{
				float result = metric.compare(BC2, b.second);
				if(result > 0.85) // 0.875... for 1 mismatch 8 bp. The only accepted mismatch rate
				{
					if(result > max2) 
					{
						multiple2 = false;
						same2 = "Barcodes with same similarity value:";
						bestMatch2 = b;
						max2 = result;
					} 
					else if(result == max2) 
					{
						multiple2 = true;
						same2 += "\t" + b.second;
					}
				}
			}
		
			if(bestMatch1 == null && bestMatch2 == null) return null;
			if(bestMatch1 == bestMatch2) 
			{
				bestMatch1.count++;
				bestMatch2.count++;
				Logger.write("[COUNTED]");
				whichReturn = 1;
			}
			else if(bestMatch1 == null && bestMatch2 != null && !multiple2)
			{
				bestMatch2.count++;
				Logger.write("[COUNTED]");
				whichReturn = 2;
			}
			else if(bestMatch2 == null && bestMatch1 != null && !multiple1)
			{
				bestMatch1.count++;
				Logger.write("[COUNTED]");
				whichReturn = 1;
			}
//			else
//			{
//				if(multiple1 && multiple2)
//				{
//					Logger.log.write("[NOT COUNTED]");
//					// They cannot be saved by the other... Stop
//				}
//				else if(multiple1)
//				{
//					// Assume barcode2 is correct, check if barcode1 can fit one of the stuffs
//					String[] possibleSave = same1.split("\t");
//					boolean flag = false;
//					for(String save:possibleSave)
//					{
//						if(bestMatch2.first.equals(save))
//						{
//							bestMatch2.count++;
//							flag = true;
//							Logger.log.write("[COUNTED]");
//							whichReturn = 2;
//							break;
//						}
//					}
//					if(!flag) Logger.log.write("[NOT COUNTED]");
//				}
//				else if(!same2.equals("Barcodes with same similarity value:"))
//				{
//					// Assume barcode1 is correct, check if barcode2 can fit one of the stuffs
//					String[] possibleSave = same2.split("\t");
//					boolean flag = false;
//					for(String save:possibleSave)
//					{
//						if(bestMatch1.second.equals(save))
//						{
//							bestMatch1.count++;
//							flag = true;
//							Logger.log.write("[COUNTED]");
//							whichReturn = 1;
//							break;
//						}
//					}
//					if(!flag) Logger.log.write("[NOT COUNTED]");
//				}
//				else
//				{
//					Logger.log.write("[NOT COUNTED]");
//					// None has other options... Stop
//				}
//			}
			
			Logger.write(BC1 + " First  barcode found.");
			if(bestMatch1 != null) 
			{
				Logger.write(" Best match = " + bestMatch1.first + " (" + bestMatch1.name + " " + bestMatch1.id + ") : " + max1);
				if(!same1.equals("Barcodes with same similarity value:")) Logger.write("[" + same1 + "]");
			}
			else Logger.write(" No Match.");
			Logger.write(", " + BC2 + " Second barcode found.");
			if(bestMatch2 != null) 
			{
				Logger.write(" Best match = " + bestMatch2.second + " (" + bestMatch2.name + " " + bestMatch2.id + ") : " + max2);
				if(!same2.equals("Barcodes with same similarity value:")) Logger.write("[" + same2 + "]");
			}
			Logger.write("\n");
			
			Parameters.overlapBoth++;
			
			switch(whichReturn)
			{
				case 1: return bestMatch1;
				case 2: return bestMatch2;
				default: return null;
			}
		}
		else if(BC1 != null) 
		{
			boolean returnBestMatch = false;
			Levenshtein metric = new Levenshtein();
			Barcode bestMatch = null;
			float max = 0;
			boolean multiple = false;
			String same = "Barcodes with same similarity value:";
			for(Barcode b:Parameters.bc)
			{
				float result = metric.compare(BC1, b.first);
				if(result > 0.9)  // 0.909... for 1 mismatch 11 bp. The only accepted mismatch rate
				{
					if(result > max) 
					{
						multiple = false;
						same = "Barcodes with same similarity value:";
						bestMatch = b;
						max = result;
					} 
					else if(result == max) 
					{
						multiple = true;
						same += "\t" + b.first;
					}
				}
			}
			
			if(!multiple && bestMatch != null)
			{
				Logger.write("[COUNTED]");
				bestMatch.count++;
				returnBestMatch = true;
			}
			else Logger.write("[NOT COUNTED]");
			
			Logger.write(BC1 + " First  barcode found.");
			if(bestMatch != null) 
			{
				Logger.write(" Best match = " + bestMatch.first + " (" + bestMatch.name + " " + bestMatch.id + ") : " + max);
				if(!same.equals("Barcodes with same similarity value:")) Logger.write("[" + same + "]");
			}
			Logger.write("\n");
			Parameters.overlapBC1++;
			
			if(returnBestMatch) return bestMatch;
			return null;
		}
		else if(BC2 != null) 
		{
			boolean returnBestMatch = false;
			Levenshtein metric = new Levenshtein();
			Barcode bestMatch = null;
			float max = 0;
			boolean multiple = false;
			String same = "Barcodes with same similarity value:";
			for(Barcode b:Parameters.bc)
			{
				float result = metric.compare(BC2, b.second);
				if(result > 0.85) // 0.875... for 1 mismatch 8 bp. The only accepted mismatch rate
				{
					if(result > max) 
					{
						multiple = false;
						same = "Barcodes with same similarity value:";
						bestMatch = b;
						max = result;
					} 
					else if(result == max) 
					{
						multiple = true;
						same += "\t" + b.second;
					}
				}
			}
			
			if(!multiple && bestMatch != null)
			{
				Logger.write("[COUNTED]");
				bestMatch.count++;
				returnBestMatch = true;
			}
			else Logger.write("[NOT COUNTED]");
			
			Logger.write(BC2 + " Second  barcode found.");
			if(bestMatch != null) 
			{
				Logger.write(" Best match = " + bestMatch.second + " (" + bestMatch.name + " " + bestMatch.id + ") : " + max);
				if(!same.equals("Barcodes with same similarity value:")) Logger.write("[" + same + "]");
			}
	
			Logger.write("\n");
			Parameters.overlapBC2++;
			
			if(returnBestMatch) return bestMatch;
			return null;
		}
		
		return null;
	}
	
	private static String getAlignedStringAtPos(SAMRecord samRecord, int start, int end)
	{
		String res = samRecord.getReadString();
		List<CigarElement> l = samRecord.getCigar().getCigarElements();
		int readStart = samRecord.getAlignmentStart();
		int readEnd = readStart; // Quality check
		boolean firstOcc = true;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					readEnd += cigar.getLength();
					break;
				case N:
					res = res.substring(0, readEnd - readStart) + PolySeq.getPolyDel(cigar.getLength()) + res.substring(readEnd - readStart, res.length());
					readEnd += cigar.getLength();
					break;
				case D:
					res = res.substring(0, readEnd - readStart) + PolySeq.getPolyDel(cigar.getLength()) + res.substring(readEnd - readStart, res.length());
					readEnd += cigar.getLength();
					break;
				case EQ:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case H:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case I:
					res = res.substring(0, readEnd - readStart) + res.substring(readEnd - readStart + 1, res.length()); // For an insertion, I "just" remove the part that is inserted
					break;
				case P:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
				case S:
					if(firstOcc) readStart -= cigar.getLength(); // S is at the beginning of the read, i.e. the beginning of the read string starts before
					else readEnd += cigar.getLength();
					break;
				case X:
					System.err.println("CIGAR = " + samRecord.getCigar());
					System.exit(-1);
			}
			firstOcc = false;
		}
		readEnd--; // Because the last letter is at the index before
		
		// Substring the read to the desired location (if overlapping, if not return null)
		if(readEnd >= end && readStart <= end) 
		{
			if(readStart >= start) // Only the beginning of the read overlaps the barcode
			{
				int i1 = 0;
				int i2 = end - readStart + 1;
				return PolySeq.POLY_N[(end - start + 1) - (i2 - i1)] + res.substring(i1, i2); // Complete with polyN if missing part of barcode
			}
			else // Barcode in full within the read
			{
				int i1 = start - readStart;
				int i2 = end - readStart + 1;
				return res.substring(i1, i2);
			}
		}
		if(readEnd <= end)
		{
			if(readStart >= start) // Read within barcode (should not be possible, unless one read is smaller than the barcode)
			{
				int i1 = 0;
				int i2 = readEnd - readStart + 1;
				return res.substring(i1, i2);
			}
			if(readEnd >= start) // Only the beginning of barcode found at the end of the read
			{
				int i1 = start - readStart;
				int i2 = readEnd - readStart + 1;
				return res.substring(i1, i2) + PolySeq.POLY_N[(end - start + 1) - (i2 - i1)];
			}
		}
		return null;
	}
}
