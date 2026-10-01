package com.hartwig.hmftools.viridian.detection.align;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Supporting information for the alignment of reads to all virus genomes.
public record AllContigsReadAlignmentsMetrics(
        double meanReadLength,
        // Reads with an alignment dropped for clipping over their contig's start or end.
        Map<ViralContig, Integer> originClippedReads,
        // Distinct reads with at least one alignment to any contig of the oncology group.
        Map<OncologyGroup, Integer> readCountsByOncologyGroup,
        // Per contig and per read, how many alignments the read had on the contig.
        Map<ViralContig, List<Integer>> alignmentCountsByContig
)
{
    public AllContigsReadAlignmentsMetrics
    {
        if(meanReadLength <= 0)
        {
            throw new IllegalArgumentException("Invalid meanReadLength: " + meanReadLength);
        }
    }
}
