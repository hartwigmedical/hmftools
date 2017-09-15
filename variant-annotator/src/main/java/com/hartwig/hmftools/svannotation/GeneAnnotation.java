package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class GeneAnnotation {

    private final Breakend breakend;
    private final String geneName;
    private final String stableId;
    private final String entrezId;
    private final int strand;
    private final List<Transcript> transcripts = Lists.newArrayList();
    private final List<String> synonyms;

    GeneAnnotation(final Breakend breakend, final String geneName, final List<String> synonyms, final String stableId,
            final String entrezId, final int strand) {
        this.breakend = breakend;
        this.geneName = geneName;
        this.synonyms = synonyms;
        this.stableId = stableId;
        this.entrezId = entrezId;
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

    public String getEntrezId() {
        return entrezId;
    }

    public int getStrand() {
        return strand;
    }

    void addTranscriptAnnotation(final Transcript a) {
        transcripts.add(a);
    }

    public List<Transcript> getTranscripts() {
        return ImmutableList.copyOf(transcripts);
    }

    public Transcript getCanonical() {
        return transcripts.stream().filter(Transcript::isCanonical).findFirst().orElse(null);
    }

    public List<String> getSynonyms() {
        return ImmutableList.copyOf(synonyms);
    }
}
