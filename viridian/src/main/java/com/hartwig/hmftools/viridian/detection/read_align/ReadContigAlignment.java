package com.hartwig.hmftools.viridian.detection.read_align;

import java.util.List;

// One read's best alignment on one contig, and how many alignments it has on that contig.
public record ReadContigAlignment(
        ViralReadAlignment best,
        int alignmentCount
)
{
    static ReadContigAlignment from(List<ViralReadAlignment> contigAlignments)
    {
        return new ReadContigAlignment(
                contigAlignments.stream().min(ViralReadAlignment.BEST_FIT_FIRST).orElseThrow(), contigAlignments.size());
    }

    public int divergence()
    {
        return best.divergence();
    }
}
