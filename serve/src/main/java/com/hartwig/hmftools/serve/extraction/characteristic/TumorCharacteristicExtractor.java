package com.hartwig.hmftools.serve.extraction.characteristic;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorCharacteristicExtractor {

    private static final Logger LOGGER = LogManager.getLogger(TumorCharacteristicExtractor.class);

    @NotNull
    private final Set<String> microsatelliteUnstableEvents;
    @NotNull
    private final Set<String> microsatelliteStableEvents;
    @NotNull
    private final Set<String> highTumorMutationalLoadEvents;
    @NotNull
    private final Set<String> lowTumorMutationalLoadEvents;
    @NotNull
    private final Set<String> highTumorMutationalBurdenEvents;
    @NotNull
    private final Set<String> lowTumorMutationalBurdenEvents;
    @NotNull
    private final Set<String> hrDeficiencyEvents;
    @NotNull
    private final Set<String> hpvPositiveEvents;
    @NotNull
    private final Set<String> ebvPositiveEvents;

    public TumorCharacteristicExtractor(@NotNull final Set<String> microsatelliteUnstableEvents,
            @NotNull final Set<String> microsatelliteStableEvents, @NotNull final Set<String> highTumorMutationalLoadEvents,
            @NotNull final Set<String> lowTumorMutationalLoadEvents, @NotNull final Set<String> highTumorMutationalBurdenEvents,
            @NotNull final Set<String> lowTumorMutationalBurdenEvents, @NotNull final Set<String> hrDeficiencyEvents,
            @NotNull final Set<String> hpvPositiveEvents, @NotNull final Set<String> ebvPositiveEvents) {
        this.microsatelliteUnstableEvents = microsatelliteUnstableEvents;
        this.microsatelliteStableEvents = microsatelliteStableEvents;
        this.highTumorMutationalLoadEvents = highTumorMutationalLoadEvents;
        this.lowTumorMutationalLoadEvents = lowTumorMutationalLoadEvents;
        this.highTumorMutationalBurdenEvents = highTumorMutationalBurdenEvents;
        this.lowTumorMutationalBurdenEvents = lowTumorMutationalBurdenEvents;
        this.hrDeficiencyEvents = hrDeficiencyEvents;
        this.hpvPositiveEvents = hpvPositiveEvents;
        this.ebvPositiveEvents = ebvPositiveEvents;
    }

    @Nullable
    public TumorCharacteristic extract(@NotNull EventType type, @NotNull String event, @NotNull String cutoff) {
        if (type == EventType.CHARACTERISTIC) {
            TumorCharacteristicAnnotation characteristic = determineCharacteristic(event);
            if (characteristic == null) {
                LOGGER.warn("Could not extract characteristic from '{}'", event);
            } else {
                TumorCharacteristicsComparator comparator = determineComparator(characteristic, cutoff);
                Double interpretedCutoff = determineCutoff(characteristic, cutoff);
                return ImmutableTumorCharacteristic.builder()
                        .annotation(characteristic)
                        .comparator(comparator)
                        .cutoff(interpretedCutoff)
                        .build();
            }
        }
        return null;
    }

    @Nullable
    @VisibleForTesting
    public TumorCharacteristicsComparator determineComparator(@Nullable TumorCharacteristicAnnotation characteristic,
            @NotNull String cutoff) {
        if (!cutoff.equals(Strings.EMPTY)) {
            String[] cutoffSplit = cutoff.split(" ");

            if (cutoff.equals("MSI high")) {
                return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
            } else if (cutoff.equals("HRD pos")) {
                return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
            } else if (cutoffSplit.length == 3) {
                if (cutoffSplit[1].equals(">=")) {
                    return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
                } else if (cutoffSplit[1].equals("<=")) {
                    return TumorCharacteristicsComparator.EQUAL_OR_LOWER;
                } else if (cutoffSplit[1].equals("<")) {
                    return TumorCharacteristicsComparator.LOWER;
                } else if (cutoffSplit[1].equals(">")) {
                    return TumorCharacteristicsComparator.GREATER;
                } else {
                    LOGGER.warn("Could not determine greater of smaller cut-off");
                    return null;
                }
            } else {
                LOGGER.warn("cutoff value '{}' couldn't be determined", cutoff);
                return null;
            }
        } else {
            if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE) {
                return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
            } else if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_STABLE) {
                return TumorCharacteristicsComparator.LOWER;
            } else if (characteristic == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD) {
                return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
            } else if (characteristic == TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD) {
                return TumorCharacteristicsComparator.LOWER;
            } else if (characteristic == TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT) {
                return TumorCharacteristicsComparator.EQUAL_OR_GREATER;
            } else {
                return null;
            }
        }
    }

    @Nullable
    @VisibleForTesting
    public Double determineCutoff(@Nullable TumorCharacteristicAnnotation characteristic, @NotNull String cutoff) {
        if (!cutoff.equals(Strings.EMPTY)) {
            String[] cutoffSplit = cutoff.split(" ");
            if (cutoff.equals("MSI high")) {
                return (double) 4;
            } else if (cutoff.equals("HRD pos")) {
                return 0.5;
            } else if (cutoffSplit.length == 3) {
                return Double.valueOf(cutoffSplit[2]);
            } else {
                LOGGER.warn("cutoff value '{}' couldn't be determined", cutoff);
                return null;
            }
        } else {
            //HMF definitions cut-off
            if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE) {
                return (double) 4;
            } else if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_STABLE) {
                return (double) 4;
            } else if (characteristic == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD) {
                return (double) 140;
            } else if (characteristic == TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD) {
                return (double) 140;
            } else if (characteristic == TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT) {
                return 0.5;
            } else {
                return null;
            }
        }
    }

    @Nullable
    private TumorCharacteristicAnnotation determineCharacteristic(@NotNull String event) {
        if (microsatelliteUnstableEvents.contains(event)) {
            return TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE;
        } else if (microsatelliteStableEvents.contains(event)) {
            return TumorCharacteristicAnnotation.MICROSATELLITE_STABLE;
        } else if (highTumorMutationalLoadEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD;
        } else if (lowTumorMutationalLoadEvents.contains(event)) {
            return TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD;
        } else if (highTumorMutationalBurdenEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN;
        } else if (lowTumorMutationalBurdenEvents.contains(event)) {
            return TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN;
        } else if (hrDeficiencyEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT;
        } else if (hpvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HPV_POSITIVE;
        } else if (ebvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.EBV_POSITIVE;
        }

        return null;
    }
}