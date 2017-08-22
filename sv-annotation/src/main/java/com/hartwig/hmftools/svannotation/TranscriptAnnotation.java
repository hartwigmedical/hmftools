package com.hartwig.hmftools.svannotation;

public class TranscriptAnnotation {

    private final GeneAnnotation gene;
    private final String transcriptId;
    private final int exonUpstream;
    private final int exonDownstream;
    private final int exonMax;
    private final boolean canonical;

    TranscriptAnnotation(final GeneAnnotation gene, final String transcriptId, final int exonUpstream, final int exonDownstream,
            final int exonMax, final boolean canonical) {
        this.gene = gene;
        this.transcriptId = transcriptId;
        this.exonUpstream = exonUpstream;
        this.exonDownstream = exonDownstream;
        this.exonMax = exonMax;
        this.canonical = canonical;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public boolean isExonic() {
        return exonUpstream > 0 && exonUpstream == exonDownstream;
    }

    public boolean isPromoter() {
        return exonUpstream == 0 && exonDownstream == 1;
    }

    public boolean isIntronic() {
        return exonUpstream > 0 && (exonDownstream - exonUpstream) == 1;
    }

    public boolean isCanonical() {
        return canonical;
    }

    public GeneAnnotation getGene() {
        return gene;
    }

    public int getExonUpstream() {
        return exonUpstream;
    }

    public int getExonDownstream() {
        return exonDownstream;
    }

    public int getExonMax() {
        return exonMax;
    }
}
