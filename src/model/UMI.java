package model;

import java.util.HashSet;

import tools.Utils;

public class UMI 
{
	public HashSet<String> umis = null;
	public int correctedSize = -1;
	
	public UMI() 
	{
		umis = new HashSet<String>();
	}
	
	public void addUMI(String umi)
	{
		umis.add(umi);
	}
	
	public int getCorrectedSize()
	{
		if(correctedSize == -1)
		{
			if(Parameters.hammingDistanceUMI == 0) { correctedSize = umis.size(); return correctedSize; } // If not sequencing error correction
			correctedSize = 0;

			if(umis.size() > 1) // UMI correction needed
			{
				String[] list = umis.toArray(new String[umis.size()]);
			    for (int i = 0; i < list.length; i++) 
			    {
			    	boolean flag = false;
			    	for (int j = i + 1; j < list.length; j++)
			    	{
			    		if(Utils.hammingDistance(list[i], list[j]) <= Parameters.hammingDistanceUMI)
			    		{
			    			flag = true;
			    			break;
			    		}
			    	}
			    	if(!flag) correctedSize++;
		        }
			}
			else correctedSize = umis.size();
		}
		return correctedSize; // Should be corrected for sequencing errors mismatches*/
	}
}
