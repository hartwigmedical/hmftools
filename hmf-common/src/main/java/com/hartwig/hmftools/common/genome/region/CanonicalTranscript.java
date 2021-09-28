package com.hartwig.hmftools.common.genome.region;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CanonicalTranscript extends TranscriptRegion {

    @NotNull
    String geneId();

    long geneStart();

    long geneEnd();

    int exonCount();

    int codingExons();

    long exonStart();

    long exonEnd();

    long exonBases();

    Strand strand();

    long codingStart();

    long codingEnd();

    int codingBases();
}
