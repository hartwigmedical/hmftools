package com.hartwig.hmftools.tars.common;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

// A transcript's alt contig and the genomic exons it maps back to.
public record ContigEntry(
        String contigName,
        int contigStart,
        int contigEnd,
        String geneId,
        String geneName,
        String transName,
        String chromosome,
        int strand,
        List<BaseRegion> exonSpans)
{
}
