package com.hartwig.hmftools.redux.splice.rescue;

// Supplementary record under consideration as a merge partner. Index is its position in the caller's
// supplementary list so the resolver can report which supps to drop without exposing SAMRecord.
public class RescueSupplementary
{
    public final int Index;
    public final String Chromosome;
    public final boolean ForwardStrand;
    public final int Start;        // 1-based
    public final String Cigar;
    public final int Mapq;

    public RescueSupplementary(
            final int index, final String chromosome, final boolean forwardStrand,
            final int start, final String cigar, final int mapq)
    {
        Index = index;
        Chromosome = chromosome;
        ForwardStrand = forwardStrand;
        Start = start;
        Cigar = cigar;
        Mapq = mapq;
    }
}
