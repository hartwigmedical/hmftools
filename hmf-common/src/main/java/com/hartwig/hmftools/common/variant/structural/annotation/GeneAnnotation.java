package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneAnnotation {

    @NotNull
    private final EnrichedStructuralVariant variant;
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
    @NotNull
    private final List<Integer> entrezIds;
    @NotNull
    private final String karyotypeBand;

    public GeneAnnotation(@NotNull final EnrichedStructuralVariant variant, final boolean isStart, @NotNull final String geneName,
            @NotNull final String stableId, final int strand, @NotNull final List<String> synonyms, @NotNull final List<Integer> entrezIds,
            @NotNull final String karyotypeBand) {
        this.variant = variant;
        this.isStart = isStart;
        this.geneName = geneName;
        this.stableId = stableId;
        this.strand = strand;
        this.synonyms = synonyms;
        this.entrezIds = entrezIds;
        this.karyotypeBand = karyotypeBand;
    }

    @NotNull
    public EnrichedStructuralVariant variant() {
        return variant;
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

    @NotNull
    public List<Integer> entrezIds() {
        return entrezIds;
    }

    @NotNull
    public String karyotypeBand() {
        return karyotypeBand;
    }
}
