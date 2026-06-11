package com.hartwig.hmftools.redux.splice.tailextend;

// Re-evaluates a lifted terminal region (bwa's M-tail anchor + soft-clip) against the contiguous
// reference, independent of bwa's tx-contig M/S split. bwa scored the alignment over the spliced
// transcript contig, so a locally-negative stretch (e.g. a 3M that is 2 matches + 1 mismatch, scoring
// -2) survived as M because the surrounding alignment carried it. Lifted back to the genome that
// decision no longer holds, so the boundary is recomputed here: the reclaim length is the prefix that
// maximises cumulative bwa-mem score walking from the near-exon boundary outward.
public final class BoundaryReclaim
{
    // bwa-mem defaults: match +1, mismatch -4.
    public static final int MATCH_SCORE = 1;
    public static final int MISMATCH_SCORE = -4;

    private BoundaryReclaim() {}

    // read and ref are aligned index-for-index, ordered from the near-exon boundary outward. Returns the
    // longest prefix in [0, len] reaching the highest cumulative score; 0 if no positive-scoring prefix.
    // >= keeps the longest tie, so an internal mismatch followed by genuine matches reclaims those matches
    // (a tie can only be reached by a match step, so the prefix never ends on a trailing mismatch).
    public static int maxScoringPrefix(final byte[] read, final byte[] ref)
    {
        final int len = Math.min(read.length, ref.length);
        int score = 0;
        int bestScore = 0;
        int bestLength = 0;
        for(int i = 0; i < len; ++i)
        {
            score += basesEqualIgnoreCase(read[i], ref[i]) ? MATCH_SCORE : MISMATCH_SCORE;
            if(score >= bestScore && score > 0)
            {
                bestScore = score;
                bestLength = i + 1;
            }
        }
        return bestLength;
    }

    // Highest cumulative score over any prefix [0, len] walking from the boundary outward (>= 0). This is the
    // value maxScoringPrefix reaches - use it to compare a contiguous extension against an alternative placement.
    public static int maxPrefixScore(final byte[] read, final byte[] ref)
    {
        final int len = Math.min(read.length, ref.length);
        int score = 0;
        int bestScore = 0;
        for(int i = 0; i < len; ++i)
        {
            score += basesEqualIgnoreCase(read[i], ref[i]) ? MATCH_SCORE : MISMATCH_SCORE;
            if(score > bestScore)
                bestScore = score;
        }
        return bestScore;
    }

    // Straight bwa-mem score of read vs ref over their overlapping length (no clipping). Used to score an anchor
    // against a fixed placement (e.g. bwa's far-exon position).
    public static int score(final byte[] read, final byte[] ref)
    {
        int total = 0;
        final int len = Math.min(read.length, ref.length);
        for(int i = 0; i < len; ++i)
            total += basesEqualIgnoreCase(read[i], ref[i]) ? MATCH_SCORE : MISMATCH_SCORE;
        return total;
    }

    // Reverses a window so a region whose near-exon boundary sits at its end can be walked boundary-outward.
    public static byte[] reversed(final byte[] bases)
    {
        final byte[] out = new byte[bases.length];
        for(int i = 0; i < bases.length; ++i)
            out[i] = bases[bases.length - 1 - i];
        return out;
    }

    static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
            return true;
        return (a & ~0x20) == (b & ~0x20);
    }
}
