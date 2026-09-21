package com.hartwig.hmftools.viridian.integration.seq_align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import htsjdk.samtools.CigarElement;

public record ViralSequenceAlignment(
        ViralContig contig,
        int position,
        Orientation orientation,
        // TODO: use real CIGAR type rather than String
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
        if(sequenceLength <= 0)
        {
            throw new IllegalArgumentException("Invalid sequenceLength: " + sequenceLength);
        }
    }

    public int alignedLength()
    {
        return cigarElementsFromStr(cigar).stream()
                .filter(element -> element.getOperator().isAlignment())
                .mapToInt(CigarElement::getLength)
                .sum();
    }
}
