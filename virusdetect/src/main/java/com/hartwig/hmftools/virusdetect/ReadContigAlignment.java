package com.hartwig.hmftools.virusdetect;

import java.util.List;

// One read's best alignment on one contig, and how many alignments it has on that contig.
public record ReadContigAlignment(
        ViralAlignment best,
        int alignmentCount
)
{
    static ReadContigAlignment from(List<ViralAlignment> contigAlignments)
    {
        return new ReadContigAlignment(
                contigAlignments.stream().min(ViralAlignment.BEST_FIT_FIRST).orElseThrow(), contigAlignments.size());
    }

    public int divergence()
    {
        return best.divergence();
    }
}
