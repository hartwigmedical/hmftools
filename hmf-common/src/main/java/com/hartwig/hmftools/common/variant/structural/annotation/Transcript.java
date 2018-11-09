package com.hartwig.hmftools.common.variant.structural.annotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transcript {

    @NotNull
    private final GeneAnnotation parent;
    @NotNull
    private final String transcriptId;
    private final int exonUpstream;
    private final int exonUpstreamPhase;
    private final int exonDownstream;
    private final int exonDownstreamPhase;
    private final int exonMax;

    private final long exonStart;
    private final long exonEnd;

    private final boolean canonical;

    @Nullable
    private final Long codingStart;

    @Nullable
    private final Long codingEnd;

    public static String TRANS_REGION_TYPE_PROMOTOR = "Promotor";
    public static String TRANS_REGION_TYPE_EXONIC = "Exonic";
    public static String TRANS_REGION_TYPE_INTRONIC = "Intronic";

    public Transcript(@NotNull final GeneAnnotation parent, @NotNull final String transcriptId,
            final int exonUpstream, final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase,
            final long exonStart, final long exonEnd,
            final int exonMax, final boolean canonical, @Nullable final Long codingStart, @Nullable final Long codingEnd)
    {
        this.parent = parent;
        this.transcriptId = transcriptId;
        this.exonUpstream = exonUpstream;
        this.exonUpstreamPhase = exonUpstreamPhase;
        this.exonDownstream = exonDownstream;
        this.exonDownstreamPhase = exonDownstreamPhase;
        this.exonMax = exonMax;
        this.exonStart = exonStart;
        this.exonEnd = exonEnd;
        this.canonical = canonical;
        this.codingStart = codingStart;
        this.codingEnd = codingEnd;
    }

    @NotNull
    public String transcriptId() {
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

    public final String getRegionType()
    {
        if(isIntronic())
            return TRANS_REGION_TYPE_INTRONIC;

        if(isExonic())
            return TRANS_REGION_TYPE_EXONIC;

        return TRANS_REGION_TYPE_PROMOTOR;
    }

    public boolean isCanonical() {
        return canonical;
    }

    @NotNull
    public GeneAnnotation parent() {
        return parent;
    }

    @NotNull
    public String geneName() {
        return parent.geneName();
    }

    public int exonUpstream() {
        return exonUpstream;
    }
    public int exonUpstreamPhase() { return exonUpstreamPhase; }
    public long exonStart() {
        return exonStart;
    }
    public long exonEnd() { return exonEnd; }

    public int exonDownstream() {
        return exonDownstream;
    }
    public int exonDownstreamPhase() {
        return exonDownstreamPhase;
    }

    public int exonMax() { return exonMax; }

    @Nullable
    public Long codingStart() {
        return codingStart;
    }

    @Nullable
    public Long codingEnd() {
        return codingEnd;
    }

    public final String toString()
    {
        return parent.geneName() + " " + transcriptId;
    }
}
