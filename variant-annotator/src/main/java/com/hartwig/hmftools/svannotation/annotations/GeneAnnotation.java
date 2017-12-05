package com.hartwig.hmftools.svannotation.annotations;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.Nullable;

public class GeneAnnotation {

    private final StructuralVariantAnnotation parent;
    private final boolean isStart;
    private final String geneName;
    private final String stableId;
    private final int strand;
    private final List<Transcript> transcripts = Lists.newArrayList();
    private final List<String> synonyms;

    public GeneAnnotation(final StructuralVariantAnnotation parent, final boolean isStart, final String geneName,
            final List<String> synonyms, final String stableId, final int strand) {
        this.parent = parent;
        this.isStart = isStart;
        this.geneName = geneName;
        this.synonyms = synonyms;
        this.stableId = stableId;
        this.strand = strand;
    }

    public StructuralVariantAnnotation getAnnotation() {
        return parent;
    }

    public StructuralVariant getVariant() {
        return parent.getVariant();
    }

    public boolean isStart() {
        return isStart;
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

    public void addTranscript(final Transcript a) {
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
