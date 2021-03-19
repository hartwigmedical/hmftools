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
    private final Set<String> microsatelliteUnstableEvents;
    @NotNull
    private final Set<String> highTumorMutationalLoadEvents;
    @NotNull
    private final Set<String> hrDeficiencyEvents;
    @NotNull
    private final Set<String> hpvPositiveEvents;
    @NotNull
    private final Set<String> ebvPositiveEvents;

    public TumorCharacteristicExtractor(@NotNull final Set<String> microsatelliteUnstableEvents,
            @NotNull final Set<String> highTumorMutationalLoadEvents, @NotNull final Set<String> hrDeficiencyEvents,
            @NotNull final Set<String> hpvPositiveEvents, @NotNull final Set<String> ebvPositiveEvents) {
        this.microsatelliteUnstableEvents = microsatelliteUnstableEvents;
        this.highTumorMutationalLoadEvents = highTumorMutationalLoadEvents;
        this.hrDeficiencyEvents = hrDeficiencyEvents;
        this.hpvPositiveEvents = hpvPositiveEvents;
        this.ebvPositiveEvents = ebvPositiveEvents;
    }

    @Nullable
    public TumorCharacteristic extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.CHARACTERISTIC) {
            TumorCharacteristic characteristic = determineCharacteristic(event);
            if (characteristic == null) {
                LOGGER.warn("Could not extract characteristic from '{}'", event);
            }
            return characteristic;
        }

        return null;
    }

    @Nullable
    private TumorCharacteristic determineCharacteristic(@NotNull String event) {
        if (microsatelliteUnstableEvents.contains(event)) {
            return TumorCharacteristic.MICROSATELLITE_UNSTABLE;
        } else if (highTumorMutationalLoadEvents.contains(event)) {
            return TumorCharacteristic.HIGH_TUMOR_MUTATIONAL_LOAD;
        } else if (hrDeficiencyEvents.contains(event)) {
            return TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT;
        } else if (hpvPositiveEvents.contains(event)) {
            return TumorCharacteristic.HPV_POSITIVE;
        } else if (ebvPositiveEvents.contains(event)) {
            return TumorCharacteristic.EBV_POSITIVE;
        }

        return null;
    }
}
