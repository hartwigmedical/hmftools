package com.hartwig.hmftools.common.region;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfTranscriptRegion implements TranscriptRegion {

    @NotNull
    public abstract String geneID();

    public abstract long geneStart();

    public abstract long geneEnd();

    public abstract long codingStart();

    public abstract long codingEnd();

    @NotNull
    public abstract Strand strand();

    @NotNull
    public abstract List<Integer> entrezId();

    @NotNull
    public abstract List<HmfExonRegion> exome();

    @Value.Derived
    @Nullable
    public HmfExonRegion exonByIndex(int index) {
        int effectiveIndex = index - 1;
        if (strand() == Strand.REVERSE) {
            // KODU: Assume the exome is sorted on genomic coordinates.
            effectiveIndex = exome().size() - index;
        }

        if (effectiveIndex >= 0 && effectiveIndex < exome().size()) {
            return exome().get(effectiveIndex);
        }

        return null;
    }
}
