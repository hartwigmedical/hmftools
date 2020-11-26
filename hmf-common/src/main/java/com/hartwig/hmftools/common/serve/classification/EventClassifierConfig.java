package com.hartwig.hmftools.common.serve.classification;

import java.util.Map;
import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventClassifierConfig {

    @NotNull
    public abstract EventPreprocessor proteinAnnotationExtractor();

    @NotNull
    public abstract String exonKeyword();

    @NotNull
    public abstract Set<String> exonRangeEvents();

    @NotNull
    public abstract Set<String> exonRangeKeywords();

    @NotNull
    public abstract Map<String, Set<String>> fusionPairAndExonRangesPerGene();

    @NotNull
    public abstract Set<String> genericGeneLevelKeywords();

    @NotNull
    public abstract Set<String> activatingGeneLevelKeywords();

    @NotNull
    public abstract Set<String> inactivatingGeneLevelKeywords();

    @NotNull
    public abstract Set<String> amplificationKeywords();

    @NotNull
    public abstract Set<String> amplificationKeyPhrases();

    @NotNull
    public abstract Set<String> deletionKeywords();

    @NotNull
    public abstract Set<String> deletionKeyPhrases();

    @NotNull
    public abstract Set<String> deletionKeywordsToSkip();

    @NotNull
    public abstract Set<String> exonicDelDupFusionEvents();

    @NotNull
    public abstract Set<String> fusionPairEventsToSkip();

    @NotNull
    public abstract Set<String> promiscuousFusionKeywords();

    @NotNull
    public abstract Set<String> signatureEvents();

    @NotNull
    public abstract Map<String, Set<String>> combinedEventsPerGene();

    @NotNull
    public abstract Map<String, Set<String>> complexEventsPerGene();
}
