package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public interface TranscriptAnnotation {

    @NotNull
    String gene();

    @NotNull
    String transcript();
}
