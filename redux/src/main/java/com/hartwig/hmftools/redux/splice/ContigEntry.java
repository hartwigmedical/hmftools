package com.hartwig.hmftools.redux.splice;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

public record ContigEntry(
        String contigName,
        String geneId,
        String geneName,
        String transName,
        String chromosome,
        // exon spans on the forward genomic strand, sorted ascending. Each span is 1-based, inclusive.
        List<BaseRegion> exonSpans)
{
    public int contigLength()
    {
        return exonSpans.stream().mapToInt(BaseRegion::baseLength).sum();
    }
}
