package com.hartwig.hmftools.viridian.integration.seq_align;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.viridian.reference.ViralContig;

public record ViralSequenceAlignment(
        ViralContig contig,
        int position,
        Orientation orientation,
        String cigar,
        int alignerScore,
        // NM tag. Edit distance over the aligned subsequence only.
        int alignedEditDistance,
        // Query sequence length.
        int sequenceLength
)
{
    public ViralSequenceAlignment
    {
        if(position < 1)
        {
            throw new IllegalArgumentException("Invalid position: " + position);
        }
        if(alignerScore < 0)
        {
            throw new IllegalArgumentException("Invalid alignerScore: " + alignerScore);
        }
        if(alignedEditDistance < 0)
        {
            throw new IllegalArgumentException("Invalid alignedEditDistance: " + alignedEditDistance);
        }
    }
}
