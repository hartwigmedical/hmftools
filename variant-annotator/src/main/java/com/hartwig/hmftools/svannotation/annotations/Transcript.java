package com.hartwig.hmftools.svannotation.annotations;

import java.sql.Struct;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class Transcript {

    private final GeneAnnotation geneAnnotation;
    private final String transcriptId;
    private final int exonUpstream;
    private final int exonUpstreamPhase;
    private final int exonDownstream;
    private final int exonDownstreamPhase;
    private final int exonMax;
    private final boolean canonical;

    public Transcript(final GeneAnnotation geneAnnotation, @NotNull final String transcriptId, final int exonUpstream,
            final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase, final int exonMax,
            final boolean canonical) {
        this.geneAnnotation = geneAnnotation;
        this.transcriptId = transcriptId;
        this.exonUpstream = exonUpstream;
        this.exonUpstreamPhase = exonUpstreamPhase;
        this.exonDownstream = exonDownstream;
        this.exonDownstreamPhase = exonDownstreamPhase;
        this.exonMax = exonMax;
        this.canonical = canonical;
    }

    @NotNull
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

    public GeneAnnotation getGeneAnnotation() {
        return geneAnnotation;
    }

    public StructuralVariant getVariant() { return geneAnnotation.getVariant(); }

    public String getGeneName() {
        return geneAnnotation.getGeneName();
    }

    public int getStrand() {
        return geneAnnotation.getStrand();
    }

    public int getExonUpstream() {
        return exonUpstream;
    }

    public int getExonUpstreamPhase() {
        return exonUpstreamPhase;
    }

    public int getExonDownstream() {
        return exonDownstream;
    }

    public int getExonDownstreamPhase() {
        return exonDownstreamPhase;
    }

    public int getExonMax() {
        return exonMax;
    }
}
