package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class GeneAnnotation {

    private final BreakendAnnotations breakend;
    private final String geneName;
    private final String stableId;
    private final int strand;
    private final List<TranscriptAnnotation> transcriptAnnotations = Lists.newArrayList();

    GeneAnnotation(final BreakendAnnotations breakend, final String geneName, final String stableId, final int strand) {
        this.breakend = breakend;
        this.geneName = geneName;
        this.stableId = stableId;
        this.strand = strand;
    }

    public BreakendAnnotations getBreakend() {
        return breakend;
    }

    public BreakendAnnotations getOtherBreakend() {
        final BreakendAnnotations start = breakend.getStructuralVariant().getStartAnnotations();
        return breakend == start ? breakend.getStructuralVariant().getEndAnnotations() : start;
    }

    public String getGeneName() {
        return geneName;
    }

    public String getStableId() {
        return stableId;
    }

    public int getStrand() {
        return strand;
    }

    void addTranscriptAnnotation(final TranscriptAnnotation a) {
        transcriptAnnotations.add(a);
    }

    public List<TranscriptAnnotation> getTranscripts() {
        return ImmutableList.copyOf(transcriptAnnotations);
    }

    public TranscriptAnnotation getCanonical() {
        return transcriptAnnotations.stream().filter(TranscriptAnnotation::isCanonical).findFirst().orElse(null);
    }
}
