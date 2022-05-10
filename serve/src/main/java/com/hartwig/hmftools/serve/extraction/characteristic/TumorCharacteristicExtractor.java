package com.hartwig.hmftools.serve.extraction.characteristic;

import java.util.Set;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorCharacteristicExtractor {

    private static final Logger LOGGER = LogManager.getLogger(TumorCharacteristicExtractor.class);

    @NotNull
    private final Set<String> microsatelliteUnstableKeyPhrases;
    @NotNull
    private final Set<String> microsatelliteStableKeyPhrases;
    @NotNull
    private final Set<String> highTumorMutationalLoadKeyPhrases;
    @NotNull
    private final Set<String> lowTumorMutationalLoadKeyPhrases;
    @NotNull
    private final Set<String> highTumorMutationalBurdenKeyPhrases;
    @NotNull
    private final Set<String> lowTumorMutationalBurdenKeyPhrases;
    @NotNull
    private final Set<String> hrDeficiencyKeyPhrases;
    @NotNull
    private final Set<String> hpvPositiveEvents;
    @NotNull
    private final Set<String> ebvPositiveEvents;

    public TumorCharacteristicExtractor(@NotNull final Set<String> microsatelliteUnstableKeyPhrases,
            @NotNull final Set<String> microsatelliteStableKeyPhrases, @NotNull final Set<String> highTumorMutationalLoadKeyPhrases,
            @NotNull final Set<String> lowTumorMutationalLoadKeyPhrases, @NotNull final Set<String> highTumorMutationalBurdenKeyPhrases,
            @NotNull final Set<String> lowTumorMutationalBurdenKeyPhrases, @NotNull final Set<String> hrDeficiencyKeyPhrases,
            @NotNull final Set<String> hpvPositiveEvents, @NotNull final Set<String> ebvPositiveEvents) {
        this.microsatelliteUnstableKeyPhrases = microsatelliteUnstableKeyPhrases;
        this.microsatelliteStableKeyPhrases = microsatelliteStableKeyPhrases;
        this.highTumorMutationalLoadKeyPhrases = highTumorMutationalLoadKeyPhrases;
        this.lowTumorMutationalLoadKeyPhrases = lowTumorMutationalLoadKeyPhrases;
        this.highTumorMutationalBurdenKeyPhrases = highTumorMutationalBurdenKeyPhrases;
        this.lowTumorMutationalBurdenKeyPhrases = lowTumorMutationalBurdenKeyPhrases;
        this.hrDeficiencyKeyPhrases = hrDeficiencyKeyPhrases;
        this.hpvPositiveEvents = hpvPositiveEvents;
        this.ebvPositiveEvents = ebvPositiveEvents;
    }

    @Nullable
    public TumorCharacteristic extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.CHARACTERISTIC) {
            TumorCharacteristicAnnotation characteristic = determineCharacteristic(event);
            if (characteristic == null) {
                LOGGER.warn("Could not extract characteristic annotation from '{}'", event);
                return null;
            }

            TumorCharacteristicsComparator comparator = determineComparator(event);
            Double interpretedCutoff = determineCutoff(comparator, event);
            return ImmutableTumorCharacteristic.builder().name(characteristic).comparator(comparator).cutoff(interpretedCutoff).build();
        }
        return null;
    }

    @Nullable
    private TumorCharacteristicAnnotation determineCharacteristic(@NotNull String event) {
        if (hasKeyPhraseMatch(event, microsatelliteUnstableKeyPhrases)) {
            return TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE;
        } else if (hasKeyPhraseMatch(event, microsatelliteStableKeyPhrases)) {
            return TumorCharacteristicAnnotation.MICROSATELLITE_STABLE;
        } else if (hasKeyPhraseMatch(event, highTumorMutationalLoadKeyPhrases)) {
            return TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD;
        } else if (hasKeyPhraseMatch(event, lowTumorMutationalLoadKeyPhrases)) {
            return TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD;
        } else if (hasKeyPhraseMatch(event, highTumorMutationalBurdenKeyPhrases)) {
            return TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN;
        } else if (hasKeyPhraseMatch(event, lowTumorMutationalBurdenKeyPhrases)) {
            return TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN;
        } else if (hasKeyPhraseMatch(event, hrDeficiencyKeyPhrases)) {
            return TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT;
        } else if (hpvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HPV_POSITIVE;
        } else if (ebvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.EBV_POSITIVE;
        }

        return null;
    }

    private static boolean hasKeyPhraseMatch(@NotNull String event, @NotNull Iterable<String> keyPhrases) {
        for (String keyPhrase : keyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }
        return false;
    }

    @Nullable
    private static TumorCharacteristicsComparator determineComparator(@NotNull String event) {
        for (TumorCharacteristicsComparator comparator : TumorCharacteristicsComparator.values()) {
            if (event.contains(comparator.keyPhrase())) {
                return comparator;
            }
        }

        return null;
    }

    @Nullable
    private static Double determineCutoff(@Nullable TumorCharacteristicsComparator comparator, @NotNull String event) {
        if (comparator == null) {
            return null;
        }

        int start = event.indexOf(comparator.keyPhrase()) + comparator.keyPhrase().length();
        return Double.parseDouble(event.substring(start).trim());
    }
}