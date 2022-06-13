package com.hartwig.hmftools.serve.extraction.immuno;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ImmunoHLAExtractor {

    public ImmunoHLAExtractor() {
    }

    private static final Logger LOGGER = LogManager.getLogger(ImmunoHLAExtractor.class);

    @Nullable
    public ImmunoHLA extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.IMMUNO_HLA) {
            String onlyHLAType = event.contains("-") ? event.split("-")[1] : event;
            String MainHlaType = onlyHLAType.contains(":") ? onlyHLAType.split(":")[0] : onlyHLAType;

            if (MainHlaType.length() != 4 ){
                LOGGER.warn("Unknown HLA type {} in knowledgebase. Investigate in more depth", MainHlaType);
            }
            return ImmutableImmunoHLA.builder().immunoHLA(event).build();
        }
        return null;
    }
}