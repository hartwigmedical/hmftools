package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_ALIGN_SCORE_MIN;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import htsjdk.samtools.Cigar;

// Alignment of an integration candidate SV's insert sequence onto a virus genome.
public record ViralInsertAlignment(
        VirusGenome genome,
        int position,
        Orientation orientation,
        Cigar cigar,
        int alignerScore,
        // NM tag. Edit distance over the aligned subsequence only.
        int alignedEditDistance,
        // Query sequence length.
        int sequenceLength
)
{
    public ViralInsertAlignment
    {
        if(cigar.isEmpty())
        {
            throw new IllegalArgumentException("Empty cigar");
        }
        if(position < 1)
        {
            throw new IllegalArgumentException("Invalid position: " + position);
        }
        if(alignerScore < INTEGRATION_ALIGN_SCORE_MIN)
        {
            throw new IllegalArgumentException("alignerScore below integration minimum: " + alignerScore);
        }
        if(alignedEditDistance < 0)
        {
            throw new IllegalArgumentException("Invalid alignedEditDistance: " + alignedEditDistance);
        }
        if(sequenceLength <= 0)
        {
            throw new IllegalArgumentException("Invalid sequenceLength: " + sequenceLength);
        }
    }

    public int alignedLength()
    {
        return sequenceLength - leftClipLength(cigar) - rightClipLength(cigar);
    }
}
