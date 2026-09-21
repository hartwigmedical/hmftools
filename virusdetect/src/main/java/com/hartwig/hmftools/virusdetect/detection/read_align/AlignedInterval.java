package com.hartwig.hmftools.virusdetect.detection.read_align;

// A contiguous run of reference bases covered by an alignment (a CIGAR M/=/X block).
public record AlignedInterval(
        // 1-based inclusive.
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
