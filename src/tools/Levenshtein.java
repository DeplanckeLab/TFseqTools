package tools;

import model.ErrorMessage;

public final class Levenshtein
{
	private final float maxCost;
	private final float insertDelete;
	private final float substitute;

	/**
	 * Constructs a new weighted Levenshtein metric. When the cost for
	 * substitution is zero Levenshtein does not satisfy the coincidence
	 * property.
	 * 
	 * @param insertDelete
	 *            positive non-zero cost of an insert or deletion operation
	 * @param substitute
	 *            positive cost of a substitute operation
	 *            
	 * @author adapted from simmetrics package
	 */
	public Levenshtein(float insertDelete, float substitute) 
	{
		if(insertDelete <= 0) new ErrorMessage("Levenshtein: insertDelete argument should be > 0");
		if(substitute < 0) new ErrorMessage("Levenshtein: substitute argument should be >= 0");
		this.maxCost = Math.max(insertDelete, substitute);
		this.insertDelete = insertDelete;
		this.substitute = substitute;
	}

	public Levenshtein() 
	{
		this(1.0f, 1.0f);
	}

	public float compare(final String a, final String b) 
	{
		if (a.isEmpty() && b.isEmpty()) return 1.0f;
		return 1.0f - (distance(a, b) / (maxCost * Math.max(a.length(), b.length())));
	}

	public float distance(final String s, final String t) 
	{
		if (s.isEmpty()) return t.length();
		if (t.isEmpty()) return s.length();
		if (s.equals(t)) return 0;

		final int tLength = t.length();
		final int sLength = s.length();

		float[] swap;
		float[] v0 = new float[tLength + 1];
		float[] v1 = new float[tLength + 1];

		// initialize v0 (the previous row of distances)
		// this row is A[0][i]: edit distance for an empty s
		// the distance is just the number of characters to delete from t
		for (int i = 0; i < v0.length; i++) v0[i] = i * insertDelete;

		for (int i = 0; i < sLength; i++) 
		{
			// first element of v1 is A[i+1][0]
			// edit distance is delete (i+1) chars from s to match empty t
			v1[0] = (i + 1) * insertDelete;

			for (int j = 0; j < tLength; j++) 
			{
				v1[j + 1] = Utils.min(v1[j] + insertDelete, v0[j + 1] + insertDelete, v0[j] + (s.charAt(i) == t.charAt(j) ? 0.0f : substitute));
			}

			swap = v0;
			v0 = v1;
			v1 = swap;
		}

		// latest results was in v1 which was swapped with v0
		return v0[tLength];
	}

	@Override
	public String toString() 
	{
		return "Levenshtein [insertDelete=" + insertDelete + ", substitute=" + substitute + "]";
	}

}