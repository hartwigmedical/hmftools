package com.hartwig.hmftools.common.genome.region;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CanonicalTranscript extends TranscriptRegion {

    @NotNull
    String geneId();

    int geneStart();

    int geneEnd();

    int exonCount();

    int codingExons();

    int exonStart();

    int exonEnd();

    int exonBases();

    Strand strand();

    int codingStart();

    int codingEnd();

    int codingBases();
}
