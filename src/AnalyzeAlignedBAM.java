import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
	public static Levenshtein metric = new Levenshtein();
	
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
					// Check duplicated read names
					String readName = samRecord.getReadName();
					if(result.containsKey(readName)) new ErrorMessage("Duplicated read names: " + readName);
					
					// Search for consistently overlapping barcode
					Barcode bc = overlapsBC(samRecord);
					
					// Add to result if found
					if(bc != null) result.put(readName, bc);
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
	
	public static ArrayList<Barcode> getBestMatchingBarcodes(String barcode, boolean isFirst)
	{
		ArrayList<Barcode> bestMatching = new ArrayList<Barcode>();
		if(barcode == null) return bestMatching;
		float max = 0;
		for(Barcode b:Parameters.bc)
		{
			String comparedTo = b.second;
			if(isFirst) comparedTo = b.first;
			float result = metric.compare(barcode, comparedTo);
			if(result > 0.9) // 0.909... for 1 mismatch 11 bp. The only accepted mismatch rate
			{
				if(result > max) 
				{
					bestMatching.clear();
					bestMatching.add(b);
					max = result;
				} 
				else if(result == max) 
				{
					bestMatching.add(b);
				}
			}
		}
		return bestMatching;
	}
	
	public static Barcode overlapsBC(SAMRecord samRecord)
	{
		// Extract barcodes
		String BC1 = getAlignedStringAtPos(samRecord, Parameters.startBC1, Parameters.endBC1); // position of barcode 1
		String BC2 = getAlignedStringAtPos(samRecord, Parameters.startBC2, Parameters.endBC2); // position of barcode 2
		if(BC1 != null && BC1.equals(Barcode.construct("-", Parameters.lBC1))) BC1 = null; // Construct of length 11
		if(BC2 != null && BC2.equals(Barcode.construct("-", Parameters.lBC2))) BC2 = null; // Construct of length 8
		if(BC1 != null && BC1.equals(Barcode.construct("N", Parameters.lBC1))) BC1 = null; // Construct of length 11
		if(BC2 != null && BC2.equals(Barcode.construct("N", Parameters.lBC2))) BC2 = null; // Construct of length 8
		
		// Get best matching barcodes
		ArrayList<Barcode> matchingBC1 = getBestMatchingBarcodes(BC1, true);
		ArrayList<Barcode> matchingBC2 = getBestMatchingBarcodes(BC2, false);
		
		// Count them or not
		// No overlap
		if(matchingBC1.isEmpty() && matchingBC2.isEmpty()) return null;
		
		// Overlap both
		if(!matchingBC1.isEmpty() && !matchingBC2.isEmpty()) 
		{
			Barcode bestMatch = null;
			
			// Intersection
			String bc1 = Barcode.toString(matchingBC1);
			matchingBC1.retainAll(matchingBC2); // Retaining barcodes (comparing object refs)
			
			// Only one in common
			if(matchingBC1.size() == 1)
			{
				bestMatch = matchingBC1.get(0);
				Logger.write("[COUNTED] Best match: BC1&BC2[" + bestMatch + "]\n");
				Parameters.overlapBoth++;
			}
			else if(matchingBC1.size() > 1)// Multiple in common
			{
				Logger.write("[NOT COUNTED] Multiple barcodes in common: BC1" + bc1 + " - BC2" + Barcode.toString(matchingBC2) + "\n");
			}
			else // No one in common
			{
				Logger.write("[NOT COUNTED] None of the found barcodes are in common: " + bc1 + " - BC2" + Barcode.toString(matchingBC2) + "\n");
			}
			return bestMatch;
		}
		
		// Overlap BC1 only
		if(!matchingBC1.isEmpty() && matchingBC2.isEmpty()) 
		{
			Barcode bestMatch = null;
			
			// Only one
			if(matchingBC1.size() == 1)
			{
				bestMatch = matchingBC1.get(0);
				Logger.write("[COUNTED] Best match: BC1[" + bestMatch + "] & BC2[Empty]\n");
				Parameters.overlapBC1++;
			}
			else if(matchingBC1.size() > 1)// Multiple
			{
				Logger.write("[NOT COUNTED] Multiple barcodes found: BC1" + Barcode.toString(matchingBC1) + " & BC2[Empty]\n");
			}
			else
			{
				// This should not happen
				new ErrorMessage("This should not happen [AnalyzeAlignedBAM, l.173]");
			}
			return bestMatch;
		}
		
		// Overlap BC2 only
		if(matchingBC1.isEmpty() && !matchingBC2.isEmpty()) 
		{
			Barcode bestMatch = null;
			
			// Only one
			if(matchingBC2.size() == 1)
			{
				bestMatch = matchingBC2.get(0);
				Logger.write("[COUNTED] Best match: BC1[Empty] & BC2[" + bestMatch + "]\n");
				Parameters.overlapBC2++;
			}
			else if(matchingBC2.size() > 1)// Multiple
			{
				Logger.write("[NOT COUNTED] Multiple barcodes found: BC1[Empty] & BC2" + Barcode.toString(matchingBC1) + "\n");
			}
			else
			{
				// This should not happen
				new ErrorMessage("This should not happen [AnalyzeAlignedBAM, l.199]");
			}
			return bestMatch;
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
