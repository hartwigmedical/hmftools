package com.hartwig.hmftools.viridian.detection.align;

import java.util.List;

// One read's best alignment on one contig, and how many alignments it has on that contig.
public record ReadContigAlignment(
        ViralReadAlignment best,
        int alignmentCount
)
{
    static ReadContigAlignment from(List<ViralReadAlignment> readContigAlignments)
    {
        // TODO: need to assert that all the alignments are for the same read and contig? or better interface to enforce this?
        return new ReadContigAlignment(
                readContigAlignments.stream().min(ViralReadAlignment.BEST_FIT_FIRST).orElseThrow(), readContigAlignments.size());
    }

    public int divergence()
    {
        return best.divergence();
    }
}
