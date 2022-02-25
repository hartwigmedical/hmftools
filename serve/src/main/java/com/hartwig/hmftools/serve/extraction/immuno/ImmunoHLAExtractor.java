package com.hartwig.hmftools.serve.extraction.immuno;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ImmunoHLAExtractor {

    public ImmunoHLAExtractor() {
    }

    @Nullable
    public ImmunoHLA extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.IMMUNO_HLA) {
            return ImmutableImmunoHLA.builder().immunoHLA(event).build();
        }
        return null;
    }
}