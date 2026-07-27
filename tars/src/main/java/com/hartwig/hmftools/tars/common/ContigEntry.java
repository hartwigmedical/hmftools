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
    // an annotation-only entry (contigStart==0, no contig): a transcript skipped from the FASTA whose exons are
    // still recorded so the liftback exon index has full ensembl coverage
    public static ContigEntry annotationOnly(
            final String geneId, final String geneName, final String transName, final String chromosome,
            final int strand, final List<BaseRegion> exonSpans)
    {
        return new ContigEntry("", 0, 0, geneId, geneName, transName, chromosome, strand, exonSpans);
    }
}
