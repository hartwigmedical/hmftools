package com.hartwig.hmftools.tars.liftback.rescue;

// Supplementary record under consideration as a merge partner. index is its position in the caller's
// supplementary list so the resolver can report which supps to drop without exposing SAMRecord.
public record RescueSupplementary(
        int index, String chromosome, boolean forwardStrand, int start, String cigar, int mapq)
{
}
