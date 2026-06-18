package com.hartwig.hmftools.tars.common;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

// One transcript segment on the alt contig. AltStart/AltEnd are 1-based inclusive positions on contigName;
// exonSpans are 1-based inclusive genomic positions. Strand (+1/-1) is forwarded to set XS:A on spliced reads.
// Annotation-only rows (altStart==0, no contig) carry the genomic exons of transcripts skipped from the FASTA
// (single-exon, all-N) so the liftback exon/junction indexes can be built from the sidecar alone; they are
// excluded from contig-coordinate lifting.
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
    public static final String ANNOTATION_ONLY_CONTIG = "annotation";

    public static ContigEntry annotationOnly(
            final String geneId, final String geneName, final String transName, final String chromosome,
            final int strand, final List<BaseRegion> exonSpans)
    {
        return new ContigEntry(ANNOTATION_ONLY_CONTIG, 0, 0, geneId, geneName, transName, chromosome, strand, exonSpans);
    }

    public boolean isAnnotationOnly()
    {
        return altStart == 0;
    }

    public int contigLength()
    {
        return altEnd - altStart + 1;
    }
}
