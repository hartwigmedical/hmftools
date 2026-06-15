package com.hartwig.hmftools.tars.common;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

// One transcript segment on the alt contig. AltStart/AltEnd are 1-based inclusive positions on contigName;
// exonSpans are 1-based inclusive genomic positions. Strand (+1/-1) is forwarded to set XS:A on spliced reads.
public record ContigEntry(
        String contigName,
        int altStart,
        int altEnd,
        String geneId,
        String geneName,
        String transName,
        String chromosome,
        int strand,
        List<BaseRegion> exonSpans)
{
    public int contigLength()
    {
        return altEnd - altStart + 1;
    }
}
