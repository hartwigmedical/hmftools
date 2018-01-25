package com.hartwig.hmftools.svannotation.annotations;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneAnnotation {

    @NotNull
    private final StructuralVariant parent;
    private final boolean isStart;
    @NotNull
    private final String geneName;
    @NotNull
    private final String stableId;
    private final int strand;
    @NotNull
    private final List<Transcript> transcripts = Lists.newArrayList();
    @NotNull
    private final List<String> synonyms;

    public GeneAnnotation(@NotNull final StructuralVariant parent, final boolean isStart, @NotNull final String geneName,
            @NotNull final String stableId, final int strand, @NotNull final List<String> synonyms) {
        this.parent = parent;
        this.isStart = isStart;
        this.geneName = geneName;
        this.stableId = stableId;
        this.strand = strand;
        this.synonyms = synonyms;
    }

    @NotNull
    public StructuralVariant variant() {
        return parent;
    }

    public boolean isStart() {
        return isStart;
    }

    @NotNull
    public String geneName() {
        return geneName;
    }

    @NotNull
    public String stableId() {
        return stableId;
    }

    public int strand() {
        return strand;
    }

    public void addTranscript(@NotNull Transcript transcript) {
        transcripts.add(transcript);
    }

    @NotNull
    public List<Transcript> transcripts() {
        return ImmutableList.copyOf(transcripts);
    }

    @Nullable
    public Transcript canonical() {
        return transcripts.stream().filter(Transcript::isCanonical).findFirst().orElse(null);
    }

    @NotNull
    public List<String> synonyms() {
        return ImmutableList.copyOf(synonyms);
    }
}
