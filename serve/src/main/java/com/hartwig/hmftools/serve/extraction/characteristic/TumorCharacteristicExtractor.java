package com.hartwig.hmftools.serve.extraction.characteristic;

import java.util.Set;

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
    private final Set<String> hrDeficiencyEvents;
    @NotNull
    private final Set<String> hpvPositiveEvents;
    @NotNull
    private final Set<String> ebvPositiveEvents;
    @NotNull
    private final Set<String> immunoHlaEvents;

    public TumorCharacteristicExtractor(@NotNull final Set<String> microsatelliteUnstableEvents,
            @NotNull final Set<String> microsatelliteStableEvents, @NotNull final Set<String> highTumorMutationalLoadEvents,
            @NotNull final Set<String> lowTumorMutationalLoadEvents, @NotNull final Set<String> hrDeficiencyEvents,
            @NotNull final Set<String> hpvPositiveEvents, @NotNull final Set<String> ebvPositiveEvents,
            @NotNull final Set<String> immunoHlaEvents) {
        this.microsatelliteUnstableEvents = microsatelliteUnstableEvents;
        this.microsatelliteStableEvents = microsatelliteStableEvents;
        this.highTumorMutationalLoadEvents = highTumorMutationalLoadEvents;
        this.lowTumorMutationalLoadEvents = lowTumorMutationalLoadEvents;
        this.hrDeficiencyEvents = hrDeficiencyEvents;
        this.hpvPositiveEvents = hpvPositiveEvents;
        this.ebvPositiveEvents = ebvPositiveEvents;
        this.immunoHlaEvents = immunoHlaEvents;
    }

    @Nullable
    public TumorCharacteristic extract(@NotNull EventType type, @NotNull String event, @NotNull String cutOff) {
        if (type == EventType.CHARACTERISTIC) {
            TumorCharacteristicAnnotation characteristic = determineCharacteristic(event);
            if (characteristic == null) {
                LOGGER.warn("Could not extract characteristic from '{}'", event);
            } else {
                String cutOffInterpretated = determineCutoff(characteristic, cutOff);
                return ImmutableTumorCharacteristic.builder().tumorCharacteristicAnnotation(characteristic).cutoff(cutOffInterpretated).build();
            }
        }
        return null;
    }

    @NotNull
    private String determineCutoff(@Nullable TumorCharacteristicAnnotation characteristic, @NotNull String cutOff) {
        if (!cutOff.equals(Strings.EMPTY)) {
            return cutOff;
        } else {
            if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE) {
                return "MSI >= 4";
            } else if (characteristic == TumorCharacteristicAnnotation.MICROSATELLITE_STABLE) {
                return "MSS < 4";
            } else if (characteristic == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD) {
                return "TML >= 140";
            } else if (characteristic == TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD) {
                return "TML < 140";
            } else if (characteristic == TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT) {
                return "HRD >= 0.5";
            } else {
                return Strings.EMPTY;
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
        } else if (hrDeficiencyEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT;
        } else if (hpvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.HPV_POSITIVE;
        } else if (ebvPositiveEvents.contains(event)) {
            return TumorCharacteristicAnnotation.EBV_POSITIVE;
        } else if (immunoHlaEvents.contains(event)) {
            return TumorCharacteristicAnnotation.IMMUNO_HLA;
        }

        return null;
    }
}
