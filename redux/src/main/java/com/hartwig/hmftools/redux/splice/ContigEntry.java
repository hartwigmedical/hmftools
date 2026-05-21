package com.hartwig.hmftools.redux.splice;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

// one transcript-shaped segment laid down on a per-chromosome alt contig. AltStart/AltEnd are 1-based inclusive
// positions on contigName; exonSpans are 1-based inclusive genomic positions on chromosome.
public record ContigEntry(
        String contigName,
        int altStart,
        int altEnd,
        String geneId,
        String geneName,
        String transName,
        String chromosome,
        List<BaseRegion> exonSpans)
{
    public int contigLength()
    {
        return altEnd - altStart + 1;
    }
}
