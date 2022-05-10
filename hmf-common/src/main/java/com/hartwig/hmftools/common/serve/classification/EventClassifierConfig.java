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
    public abstract Set<String> exonIdentifiers();

    @NotNull
    public abstract Set<String> exonKeywords();

    @NotNull
    public abstract Set<String> exonBlacklistKeyPhrases();

    @NotNull
    public abstract Set<String> specificExonEvents();

    @NotNull
    public abstract Map<String, Set<String>> fusionPairAndExonsPerGene();

    @NotNull
    public abstract Set<String> geneLevelBlacklistKeyPhrases();

    @NotNull
    public abstract Set<String> geneWildTypesKeyPhrases();

    @NotNull
    public abstract Set<String> genericGeneLevelKeyPhrases();

    @NotNull
    public abstract Set<String> activatingGeneLevelKeyPhrases();

    @NotNull
    public abstract Set<String> inactivatingGeneLevelKeyPhrases();

    @NotNull
    public abstract Set<String> amplificationKeywords();

    @NotNull
    public abstract Set<String> amplificationKeyPhrases();

    @NotNull
    public abstract Set<String> deletionBlacklistKeyPhrases();

    @NotNull
    public abstract Set<String> deletionKeywords();

    @NotNull
    public abstract Set<String> deletionKeyPhrases();

    @NotNull
    public abstract Set<String> exonicDelDupFusionKeyPhrases();

    @NotNull
    public abstract Set<String> exonicDelDupFusionEvents();

    @NotNull
    public abstract Set<String> fusionPairEventsToSkip();

    @NotNull
    public abstract Set<String> promiscuousFusionKeyPhrases();

    @NotNull
    public abstract Set<String> microsatelliteStableKeyPhrases();

    @NotNull
    public abstract Set<String> microsatelliteUnstableKeyPhrases();

    @NotNull
    public abstract Set<String> highTumorMutationalLoadKeyPhrases();

    @NotNull
    public abstract Set<String> lowTumorMutationalLoadKeyPhrases();

    @NotNull
    public abstract Set<String> highTumorMutationalBurdenKeyPhrases();

    @NotNull
    public abstract Set<String> lowTumorMutationalBurdenKeyPhrases();

    @NotNull
    public abstract Set<String> hrDeficiencyKeyPhrases();

    @NotNull
    public abstract Set<String> hlaEvents();

    @NotNull
    public abstract Set<String> hpvPositiveEvents();

    @NotNull
    public abstract Set<String> ebvPositiveEvents();

    @NotNull
    public abstract Map<String, Set<String>> combinedEventsPerGene();

    @NotNull
    public abstract Map<String, Set<String>> complexEventsPerGene();
}