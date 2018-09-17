package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.region.Strand;
import com.hartwig.hmftools.common.region.TranscriptRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CanonicalTranscript extends TranscriptRegion {

    @NotNull
    String geneID();

    long geneStart();

    long geneEnd();

    int exons();

    int codingExons();

    long exonStart();

    long exonEnd();

    long exonBases();

    Strand strand();

    long codingStart();

    long codingEnd();

    int codingBases();
}
