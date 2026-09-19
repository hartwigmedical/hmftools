package com.hartwig.hmftools.virusdetect;

// A contiguous run of reference bases covered by an alignment (a CIGAR M/=/X block).
// 1-based reference start.
public record AlignedInterval(
        int referenceStart,
        int length
)
{
    public AlignedInterval
    {
        if(referenceStart < 1)
        {
            throw new IllegalArgumentException("Invalid reference start: " + referenceStart);
        }
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid interval length: " + length);
        }
    }
}
