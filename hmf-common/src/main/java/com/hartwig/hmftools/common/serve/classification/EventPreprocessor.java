package com.hartwig.hmftools.common.serve.classification;

import org.jetbrains.annotations.NotNull;

public interface EventPreprocessor {

    @NotNull
    String apply(@NotNull String event);
}
