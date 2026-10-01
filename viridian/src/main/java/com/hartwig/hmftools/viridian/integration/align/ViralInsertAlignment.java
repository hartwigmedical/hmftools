package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_ALIGN_SCORE_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

// Alignment of a candidate integration variant's insert sequence onto a virus genome.
public record ViralInsertAlignment(
        ViralContig contig,
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
        return cigar.getCigarElements().stream()
                .filter(element -> element.getOperator().isAlignment())
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    public double scorePerAlignedBase()
    {
        int alignedLength = alignedLength();
        return alignedLength > 0 ? (double) alignerScore / alignedLength : 0.0;
    }

    public boolean passesFilters()
    {
        return alignerScore >= INTEGRATION_ALIGN_SCORE_MIN && scorePerAlignedBase() >= INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN;
    }
}
