package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class GeneAnnotation {

    private final Breakend breakend;
    private final String geneName;
    private final String stableId;
    private final int strand;
    private final List<Transcript> transcripts = Lists.newArrayList();
    private final List<String> synonyms;

    public GeneAnnotation(final Breakend breakend, final String geneName, final List<String> synonyms, final String stableId, final int strand) {
        this.breakend = breakend;
        this.geneName = geneName;
        this.synonyms = synonyms;
        this.stableId = stableId;
        this.strand = strand;
    }

    public Breakend getBreakend() {
        return breakend;
    }

    public Breakend getOtherBreakend() {
        final Breakend start = breakend.getStructuralVariant().getStart();
        return breakend == start ? breakend.getStructuralVariant().getEnd() : start;
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

    public void addTranscriptAnnotation(final Transcript a) {
        transcripts.add(a);
    }

    public List<Transcript> getTranscripts() {
        return ImmutableList.copyOf(transcripts);
    }

    @Nullable
    public Transcript getCanonical() {
        return transcripts.stream().filter(Transcript::isCanonical).findFirst().orElse(null);
    }

    public List<String> getSynonyms() {
        return ImmutableList.copyOf(synonyms);
    }
}
