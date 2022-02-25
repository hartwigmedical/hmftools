package com.hartwig.hmftools.serve.extraction.immuno;

import java.util.Set;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.immuno.ImmutableActionableHLA;
import com.hartwig.hmftools.serve.extraction.characteristic.ImmutableTumorCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsAtLeast;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ImmunoHLAExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ImmunoHLAExtractor.class);

    @NotNull
    private final Set<String> immunoHlaEvents;

    public ImmunoHLAExtractor(@NotNull final Set<String> immunoHlaEvents) {
        this.immunoHlaEvents = immunoHlaEvents;
    }

    @Nullable
    public ImmunoHLA extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.IMMUNO_HLA) {
            return ImmutableImmunoHLA.builder().immunoHLA(event).build();
        }
        return null;
    }
}
